#!/usr/bin/perl -w
# This script calculates D stats (ABBA-BABA test, or D test) and fd stats for sliding window from VCF file.
# Used for introgression analysis.
# However, this script takes up a lot of memory when it runs, which should be improved.
# And, the script does not calculate Z scores, which should also be improved.

use Bio::DB::HTS::VCF;
use Getopt::Long;
use List::Util qw(sum);
use Parallel::ForkManager;

### parse parameters
GetOptions(
	'vcf|v:s' => \$vcf_file,
	'trios|t:s' => \$trios_file,
	'map|m:s' => \$sample_pop_file,
	'interval|i:s' => \$interval_file,
	'window|w:i' => \$window,
	'step|s:i' => \$step,
	'processor|p:i' => \$processor,
	'help|h!' => \$h
);
$window ||= 30000;
$step ||= 10000;
$processor ||= 4;

$help = <<"EOF";
    Calculate D (ABBA-BABA) and fd statistics for sliding window from VCF file. The columns of input files should be tab-separated.

    Usage: perl $0 [options]
        -v   vcf file. bgzip compressed vcf with tbi format index. 
        -t   trios file. Three colums: pop1, pop2, pop3. No outgroup, outgroup is specified in map file.
        -m   map flie. Two columns: sample name, population name. Population name of outgroup sample must be 'Outgroup'.
        -i   interval file. Three columns: chromosome name, 1-based start coordinate, 1-based end coordinate.
        -w   window size for sliding window. Default is 30,000.
        -s   step size for sliding window. Default is 10,000.
        -p   number of processors to use. Default is 4.
        -h   help message.
EOF
if (!($vcf_file and $trios_file and $sample_pop_file and $interval_file) or $h) {
	print STDERR "\n$help\n";
	exit;
}

### read sample_pop_file
open SMP, '<', $sample_pop_file;
%pop_sample = ();
print STDERR "Reading map file...\n";
while (<SMP>) {
	next if /^#/;
	chomp;
	@l = split /\s+/;
	push @{$pop_sample{$l[1]}}, $l[0];
}
close SMP;

### read trios_file
print STDERR "Reading trios file...\n";
open TRIOS, '<', $trios_file;
while (<TRIOS>) {
	next if /^#/;
	chomp;
	push @trios, $_;
}
close TRIOS;

### read interval_file and generate regions for sliding window D and fd calculation.
print STDERR "Reading interval file...\n";
open INTERVAL, '<', $interval_file;
while (<INTERVAL>) {
	next if /^#/;
	chomp;
	($chr, $chr_start, $chr_end) = split /\s+/;
	$region = "$chr:$chr_start-$chr_end";
	push @regions, $region;
}
%count = ();
@regions = grep {++$count{$_} < 2} @regions;
close INTERVAL;

### create output file for every trio
for $trio (@trios) {
	@out_file = ((split /\s+/,$trio), $window, $step);
	$out_file = join '.',@out_file;
	open OUT, '>', $out_file;
	print OUT "Pop1\tPop2\tPop3\tChr\tStart\tEnd\tD\tfd\n";
}
open ALL, '>', "all.D.fd.stats";
print ALL "Pop1\tPop2\tPop3\tD\tfd\n";
close ALL;

### prepare for multi threads
$trios_len = scalar @trios;
$a = sprintf "%d", $trios_len / $processor;
$b = $trios_len % $processor;
@trios_num_for_every_processor = ();
push @trios_num_for_every_processor, $a for (1..$processor);
for (0..($b-1)) {$trios_num_for_every_processor[$_]++};
for (@trios_num_for_every_processor) {
	push @trios2, [splice @trios, 0, $_];
};

### calculate D and fd by sliding window and genome wide.
### And output D and fd of windows for every trio and genome wide.
$pm = Parallel::ForkManager->new($processor);
for $temp_trios (@trios2) {
	$pm->start and next;
	$vcf = Bio::DB::HTS::VCF->new(filename => $vcf_file);
	$header = $vcf->header();
	$sample_names = $header->get_sample_names(); #$sample_names is a ref of a array of sample names in vcf file.
	@sample_names = @$sample_names;
	@trios = @$temp_trios;
	%ABBA_total = ();
	%BABA_total = ();
	%denominator_d_total = ();
	%denominator_fd_total = ();
	for $region (@regions) {
		$iter = $vcf->query("$region");
		@region_line = split /:|-/, $region;
		$window_start = $region_line[1];
		$window_end = $region_line[1] + $window - 1;
		($A1, $A2, $A3, $A4, $B1, $B2, $B3, $B4, $t1, $t2, $t3, $t4) = (0,0,0,0,0,0,0,0,0,0,0,0);
		%ABBA = ();
		%BABA = ();
		%denominator_fd = ();
		$last_pos = $window_start - 1;
		print STDERR "Prossing region ${region}...\n";
		while ($row = $iter->next()) {
			$chr = $row->chromosome($header);
			$pos = $row->position();
			$num_to_fill = $pos - $last_pos - 1;
			$genotypes = $row->get_genotypes($header); #$genotypes is a ref of a array of integers representing genotype records.
			%sample_genotype = ();
			$ploidy = (scalar @$genotypes) / (scalar @sample_names);
			$ref = $row->reference();
			$alleles = $row->get_alleles();
			@all_genotypes = ($ref, @$alleles);
			%count = ();
			$number_genotype = scalar (grep {! $count{$_}++} @all_genotypes);
			for $i (0..$#sample_names) {
				$g1 = $i * $ploidy;
				$g2 = $g1 + $ploidy - 1;
				@{$sample_genotype{$sample_names[$i]}} = @{$genotypes}[$g1..$g2];
			}
			if ($pos >= $window_start and $pos < $window_end) {
				&process_window_line();
			}
			elsif ($pos == $window_end) {
				&process_window_line();
				&print_window_d_fd;
				&splice_d_fd_for_sliding;
				$window_start += $step;
				$window_end += $step;
			}
			elsif ($pos > $window_end) {
				$num_to_fill = $window_end - $last_pos;
				for $i (0..$#trios) {
					@temp = ();
					@temp2 = ();
					$trio = $trios[$i];
					@trio = split /\s+/, $trio;
					if ($num_to_fill) {
						for (1..$num_to_fill) {
							&push_zero();
						}
					}
				}
				&print_window_d_fd;
				&splice_d_fd_for_sliding;
				$window_start += $step;
				$window_end += $step;
				while ($pos > $window_end) {
					$window_start = $window_start + $step;
					$window_end = $window_end + $step;
					$num_to_fill = $step;
					for $i (0..$#trios) {
						@temp = ();
						@temp2 = ();
						$trio = $trios[$i];
						@trio = split /\s+/, $trio;
						if ($num_to_fill) {
							for (1..$num_to_fill) {
								&push_zero();
							}
						}
					}
					&print_window_d_fd;
					&splice_d_fd_for_sliding;
				}
				$num_to_fill = $pos - ($window_end - $step) - 1;
				&process_window_line();
			}
			$last_pos = $pos;
		}
		$key = $trios[0];
		&print_window_d_fd if ((scalar @{$ABBA{$key}}) > ($window * 0.9));
	}

	### output genome-wide D and fd
	open ALL, '>>', "all.D.fd.stats";
	for $trio (@trios) {
		$numerator_d_total = $ABBA_total{$trio} - $BABA_total{$trio};
		$denominator_d_total = $ABBA_total{$trio} + $BABA_total{$trio};
		$denominator_fd_total = $denominator_fd_total{$trio};
		$d_total = $denominator_d_total ? ($numerator_d_total / $denominator_d_total) : 0;
		$fd_total = $denominator_fd_total ? ($numerator_d_total / $denominator_fd_total) : 0;
		$trio_out = join("\t",split(/\s+/, $trio));
		print ALL "$trio_out\t$d_total\t$fd_total\n";
	}
	$pm->finish;
}
$pm->wait_all_children;
print STDERR "Output genome-wide stats to all.D.fd.stats file.\n";
print STDERR "All complete!\n";

##### subroutines #####
sub splice_d_fd_for_sliding {
	for $key (keys %ABBA) {
		print scalar @{$ABBA{$key}}; print "\n";
		print scalar @{$BABA{$key}}; print "\n";
		print scalar @{$denominator_fd{$key}}; print "\n";
		splice @{$ABBA{$key}}, 0, $step;
		splice @{$BABA{$key}}, 0, $step;
		splice @{$denominator_fd{$key}}, 0, $step;
		print scalar @{$ABBA{$key}}; print "\n";
		print scalar @{$BABA{$key}}; print "\n";
		print scalar @{$denominator_fd{$key}}; print "\n";
	}
}

sub print_window_d_fd {
	for $key (keys %ABBA) {
		@out_file = ((split /\s+/,$key), $window, $step);
		$out_file = join '.',@out_file;
		open OUT, '>>', $out_file;
		$ABBA_window = sum @{$ABBA{$key}};
		$BABA_window = sum @{$BABA{$key}};
		$numerator_d_window = $ABBA_window - $BABA_window;
		$denominator_d_window = $ABBA_window + $BABA_window;
		$denominator_fd_window = sum @{$denominator_fd{$key}};
		$d_window = $denominator_d_window ? ($numerator_d_window / $denominator_d_window) : 0;
		$fd_window = $denominator_fd_window ? ($numerator_d_window / $denominator_fd_window) : 0;
		$trio_out = join("\t",(split /\s+/,$key));
		print OUT "$trio_out\t$chr\t$window_start\t$window_end\t$d_window\t$fd_window\n";
	}
#	print STDERR "Processing for window $chr:$window_start-$window_end is finished.\n";
}

sub process_window_line {
	for $i (0..$#trios) {
		@temp = ();
		@temp2 = ();
		$trio = $trios[$i];
		@trio = split /\s+/, $trio;
		&get_trio_D_fd();
		if ($num_to_fill) {
			for (1..$num_to_fill) {
				&push_zero();
			}
		}
		&push_stats();
		&push_all();
	}
}

sub get_trio_D_fd {
	for $p (@trio) {
		for $s (@{$pop_sample{$p}}) {
			push @temp, @{$sample_genotype{$s}};
		}
	}
	for $s (@{$pop_sample{'Outgroup'}}) {
		push @temp2, @{$sample_genotype{$s}};
	}
	push @temp, @temp2;
	%count = ();
	$len1 = grep {++$count{$_} < 2} @temp;
	%count = ();
	$len2 = grep {++$count{$_} < 2} @temp2;
	if ($len1 != 2 or $len2 != 1 or $number_genotype != 2) {
		$ABBA = 0;
		$BABA = 0;
		$denominator_fd = 0;
	}
	else {
		$outgroup_name = $pop_sample{'Outgroup'}->[0];
		$outgroup_genotype = $sample_genotype{$outgroup_name}->[0];
		($A1, $A2, $A3, $A4, $B1, $B2, $B3, $B4, $t1, $t2, $t3, $t4) = (0,0,0,0,0,0,0,0,0,0,0,0);
		$AABA = 0;
		$ABAA = 0;
		for $s (@{$pop_sample{$trio[0]}}) {
			for $g (@{$sample_genotype{$s}}) {
				$B1++ if $g != $outgroup_genotype;
				$t1++;
			}
		}
		$B1 /= $t1;
		$A1 = 1 - $B1;
		for $s (@{$pop_sample{$trio[1]}}) {
			for $g (@{$sample_genotype{$s}}) {
				$B2++ if $g != $outgroup_genotype;
				$t2++;
			}
		}
		$B2 /= $t2;
		$A2 = 1 - $B2;
		for $s (@{$pop_sample{$trio[2]}}) {
			for $g (@{$sample_genotype{$s}}) {
				$B3++ if $g != $outgroup_genotype;
				$t3++;
			}
		}
		$B3 /= $t3;
		$A3 = 1 - $B3;
		for $s (@{$pop_sample{'Outgroup'}}) {
			for $g (@{$sample_genotype{$s}}) {
				$B4++ if $g != $outgroup_genotype;
				$t4++;
			}
		}
		$B4 /= $t4;
		$A4 = 1 - $B4;
		$ABBA = $A1 * $B2 * $B3 * $A4;
		$BABA = $B1 * $A2 * $B3 * $A4;
		$AABA = $A1 * $A2 * $B3 * $A4;
		$ABAA = $A1 * $B2 * $A3 * $A4;
		$ABBA1224 = $A1 * $B2 * $B2 * $A4;
		$BABA1224 = $B1 * $A2 * $B2 * $A4;
		$ABBA1334 = $A1 * $B3 * $B3 * $A4;
		$BABA1334 = $B1 * $A3 * $B3 * $A4;
		if ($t1 == 1 and $t2 == 1 and $t3 == 1 and $t4 == 1) {
			$denominator_fd = $ABBA + $AABA + $ABAA;
		}
		elsif ($B2 >= $B3) {
			$denominator_fd = $ABBA1224 - $BABA1224;
		}
		else {
			$denominator_fd = $ABBA1334 - $BABA1334;
		}
	}
}

sub push_zero {
	push @{$ABBA{$trio}}, 0;
	push @{$BABA{$trio}}, 0;
	push @{$denominator_fd{$trio}}, 0;
}

sub push_stats {
	push @{$ABBA{$trio}}, $ABBA;
	push @{$BABA{$trio}}, $BABA;
	push @{$denominator_fd{$trio}}, $denominator_fd;
}

sub push_all {
	$ABBA_total{$trio} += $ABBA;
	$BABA_total{$trio} += $BABA;
	$denominator_fd_total{$trio} += $denominator_fd;
}
