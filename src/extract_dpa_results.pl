#!/usr/bin/env perl
# Extract DPA results: residue lists and top-scoring DPA points
# Usage: perl extract_dpa_results.pl <pdbid>

#use strict;
use warnings;

my $pdbid = $ARGV[0] || die "Usage: $0 <pdbid>\n";

if (! exists $ENV{DPA_HOME}) {
  print "Must set DPA_HOME to repository directory\n";
  exit(1);
}

my $home=$ENV{DPA_HOME};

# Residue ids from dpa.pl are "chain:resSeq" tokens (e.g. "A:69") since a
# structure's chains can share residue numbering.
sub sort_chain_resid {
    return sort {
        my ($ca, $ra) = split /:/, $a, 2;
        my ($cb, $rb) = split /:/, $b, 2;
        $ca cmp $cb or $ra <=> $rb
    } @_;
}

my @dpa_site_colors = ('orange', 'yellow', 'pink');

sub pymol_resi_selection {
    my %by_chain;
    foreach my $r (@_) {
        my ($chain, $resid) = split /:/, $r, 2;
        push @{$by_chain{$chain}}, $resid;
    }
    my @clauses;
    foreach my $chain (sort keys %by_chain) {
        my @resids = sort { $a <=> $b } @{$by_chain{$chain}};
        push @clauses, "(chain $chain and resi " . join("+", @resids) . ")";
    }
    return join(" or ", @clauses);
}

# Read _ana.res3 for residue lists
print "=" x 70, "\n";
print "DPA ANALYSIS RESULTS FOR $pdbid\n";
print "=" x 70, "\n\n";

open(my $pmlf, ">", "$home/scratch/${pdbid}_pymol.pml");

if (-e "$home/structuredata/${pdbid}.pdb") {
    print $pmlf "load $home/structuredata/${pdbid}.pdb,${pdbid}_prot\n";
}

if (-e "$home/structuredata/ligand_${pdbid}.pdb") {
    print $pmlf "load $home/structuredata/ligand_${pdbid}.pdb,${pdbid}_het\n";
    print $pmlf "hide everything,${pdbid}_het\n";
    print $pmlf "show spheres, ${pdbid}_het\n";
}
if (-e "$home/sscratch/${pdbid}_topdpacluster.pdb") {
    print $pmlf "load $home/scratch/${pdbid}_topdpacluster.pdb,${pdbid}_topdpa\n";
    print $pmlf "hide everything,${pdbid}_topdpa\n";
    print $pmlf "show spheres, ${pdbid}_topdpa\n";
    print $pmlf "color magenta, ${pdbid}_topdpa\n";
    print $pmlf "toggle everything, ${pdbid}_topdpa\n";
}
if (-e "$home/scratch/${pdbid}_alldpa.pdb") {
    print $pmlf "load $home/scratch/${pdbid}_alldpa.pdb,${pdbid}_dpa\n";
    print $pmlf "hide everything,${pdbid}_dpa\n";
    print $pmlf "show spheres, ${pdbid}_dpa\n";
    print $pmlf "spectrum b, rainbow, ${pdbid}_dpa\n";
    print $pmlf "set sphere_scale, 0.3, ${pdbid}_dpa\n";
    print $pmlf "toggle everything, ${pdbid}_dpa\n";
}

if (-e "_ana.res3") {
    open(my $fh, "<", "_ana.res3") || die "Cannot open _ana.res3: $!\n";
    print "BINDING SITE RESIDUES:\n";
    print "-" x 70, "\n";

    while (<$fh>) {
        if (/PLG_BINDINGSITES>\s+(\S+)\s+---\s+(\d+)\s+--\s+(.+)$/) {
            my ($id, $site_num, $residues) = ($1, $2, $3);
            my @res_list = split /\s+/, $residues;
            print "Ligand Binding Site $site_num (from protein-ligand contacts):\n";
            print "  Residues: ", join(", ", sort_chain_resid(@res_list)), "\n";
            print "  Count: ", scalar(@res_list), "\n\n";
            my $ressel = pymol_resi_selection(@res_list);
            print $pmlf "select lig_site_${site_num}, ${pdbid}_prot and ($ressel) \n";
        }
        elsif (/DPA_BINDINGSITES>\s+(\S+)\s+--\s+(\d+)\s+--\s+(\d+)\s+--\s+(.+)$/) {
            my ($id, $site_num, $count, $residues) = ($1, $2, $3, $4);
            my @res_list = split /\s+/, $residues;
            print "DPA Predicted Site $site_num:\n";
            print "  Residues: ", join(", ", sort_chain_resid(@res_list)), "\n";
            print "  Count: $count\n\n";
            my $ressel = pymol_resi_selection(@res_list);
            print $pmlf "select dpa_site_${site_num}, ${pdbid}_prot and ($ressel) \n";
            my $color = $dpa_site_colors[($site_num - 1) % scalar(@dpa_site_colors)];
            print $pmlf "color $color, dpa_site_${site_num}\n";
            print $pmlf "show sticks, dpa_site_${site_num}\n";
        }
    }
    close($fh);
}

# Check for top DPA cluster PDB file (in scratch directory or current)
my $cluster_pdb = "${pdbid}_topdpacluster.pdb";
$cluster_pdb = "$home/scratch/$cluster_pdb" if (!-e $cluster_pdb && -e "$home/scratch/$cluster_pdb");

if (-e $cluster_pdb) {
    print "\n", "=" x 70, "\n";
    print "TOP DPA POINTS (from $cluster_pdb):\n";
    print "-" x 70, "\n";

    open(my $fh, "<", $cluster_pdb) || die "Cannot open $cluster_pdb: $!\n";
    my %clusters;
    my $point_count = 0;

    while (<$fh>) {
        # PDB format: fixed columns
        # HETATM    1  CA  TOP O   0      14.914  55.757  60.341  1.00140.09      MING   0
        if (/^HETATM/) {
            my $chain = substr($_, 21, 1);
            my $resseq = substr($_, 22, 4);
            my $x = substr($_, 30, 8);
            my $y = substr($_, 38, 8);
            my $z = substr($_, 46, 8);
            my $temp = substr($_, 60, 6);

            # Remove leading/trailing whitespace
            $resseq =~ s/^\s+|\s+$//g;
            $x =~ s/^\s+|\s+$//g;
            $y =~ s/^\s+|\s+$//g;
            $z =~ s/^\s+|\s+$//g;
            $temp =~ s/^\s+|\s+$//g;

            # Chain ID indicates cluster: O=1, P=2, Q=3, etc.
            my $cluster_id = ord($chain) - ord('O') + 1;
            $cluster_id = 0 if $chain eq 'X';  # Noise cluster

            push @{$clusters{$cluster_id}}, {
                id => $resseq,
                x => $x,
                y => $y,
                z => $z,
                score => $temp
            };
            $point_count++;
        }
    }
    close($fh);

    foreach my $cid (sort {$a <=> $b} keys %clusters) {
        my @points = @{$clusters{$cid}};
        print "\nCluster $cid: ", scalar(@points), " points\n";
        printf "  %-8s  %-10s %-10s %-10s  %-10s\n", "Point", "X", "Y", "Z", "DPA Score";

        # Sort by DPA score (highest first) and show top 10
        my @sorted = sort { $b->{score} <=> $a->{score} } @points;
        my $show = scalar(@sorted) > 10 ? 10 : scalar(@sorted);

        for (my $i = 0; $i < $show; $i++) {
            printf "  %-8d  %10.3f %10.3f %10.3f  %10.3f\n",
                $sorted[$i]{id}, $sorted[$i]{x}, $sorted[$i]{y},
                $sorted[$i]{z}, $sorted[$i]{score};
        }

        if (scalar(@sorted) > 10) {
            print "  ... (", scalar(@sorted) - 10, " more points)\n";
        }
    }

    print "\nTotal: $point_count top-scoring DPA points\n";
} else {
    print "\n", "=" x 70, "\n";
    print "NOTE: Cluster PDB file not found.\n";
    print "To generate it, run dpa.pl with the -wcpdb flag:\n";
    print "  perl dpa.pl -p $pdbid -topp 0.96 -cutoff 6 6 -wcpdb\n";
}

# Read fitting parameters
if (-e "_ana.res") {
    print "\n", "=" x 70, "\n";
    print "EXTREME VALUE DISTRIBUTION PARAMETERS:\n";
    print "-" x 70, "\n";

    open(my $fh, "<", "_ana.res") || die "Cannot open _ana.res: $!\n";
    while (<$fh>) {
        if (/TOPP\s+([\d\.]+)\s+FITTING\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)/) {
            my ($percentile, $param_a, $param_b, $std_err, $corr) = ($1, $2, $3, $4, $5);
            print "  Top percentile: $percentile\n";
            print "  EVD parameter a: $param_a\n";
            print "  EVD parameter b: $param_b\n";
            print "  Standard error: $std_err\n";
            print "  Correlation: $corr\n";
        }
    }
    close($fh);
}
close($pmlf);
print "\n", "=" x 70, "\n";
