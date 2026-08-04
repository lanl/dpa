#!/usr/bin/perl

&parse_arg();

$timetile = &timetitle();
open(prcf,"> _ana.procc");
print prcf "Record DPA analysis at $timetile\n";
print prcf "COMMAND: [";for($i=0;$i<=$#ARGV;$i++){print prcf " $ARGV[$i]";}
print prcf " -topp $toppercent -nsec $ndpasection -clu $epislon_dpa $MinPts_dpa $epislon_cluster_dpa";
print prcf "]\n\n"; 
open(rcdf,"> _ana.res");
open(rcdf2,"> _ana.res2");
open(rcdf3,"> _ana.res3");
print rcdf "RECORDING ==DPA== RESULTS at $timetile\n";
print rcdf "COMMAND: [";for($i=0;$i<=$#ARGV;$i++){print rcdf " $ARGV[$i]";}
print rcdf " -topp $toppercent -nsec $ndpasection -clu $epislon_dpa $MinPts_dpa $epislon_cluster_dpa";
print rcdf "]\n\n";
print rcdf2 "RECORDING ==DPA== RESULTS at $timetile\n";
print rcdf2 "COMMAND: [";for($i=0;$i<=$#ARGV;$i++){print rcdf2 " $ARGV[$i]";}
print rcdf2 " -topp $toppercent -nsec $ndpasection -clu $epislon_dpa $MinPts_dpa $epislon_cluster_dpa";
print rcdf2 "]\n\n";
print rcdf3 "RECORDING ==DPA== RESULTS at $timetile\n";
print rcdf3 "COMMAND: [";for($i=0;$i<=$#ARGV;$i++){print rcdf3 " $ARGV[$i]";}
print rcdf3 " -topp $toppercent -nsec $ndpasection -clu $epislon_dpa $MinPts_dpa $epislon_cluster_dpa";
print rcdf3 "]\n\n";

my @idlist=();undef @idlist;
if($pdblist eq "-p"){
    push (@idlist,$pdbid)}
else{
    open(listf,"<$pdblistf");
    $checklist=0;
    while(<listf>){             # For each PDBID list -- TOPPEST list
	$checklist++;
	next if(($checklist-$liststart)%$liststep ne 0 or $checklist<$liststart); 
	my @tmp=split /\s+/,$_; my $pdbid=@tmp[1];
	$pdbid=@tmp[0] if(@tmp[0] ne '' and @tmp[0] ne ' ');
	push (@idlist,$pdbid);
	}
    close(listf)
    }

my $Nlgclutt=$Ndpaclutt=$Npredtt=0;

foreach $pdbid ( @idlist){
    my $pdbfile=$pdbid.".pdb";my $ligandf="ligand_".$pdbid.".pdb";
#
#-- C_alpha DATA
    if(! -e "$structuredatadir/$pdbfile"){
	print prcf "$pdbid -9999 --  Not detect PROTEIN PDBfile\n";
	print "$pdbid -9999 --  Not detect PROTEIN PDBfile\n";}
    my ($nca,$cacrd)=&read_pdb($structuredatadir,$pdbfile,$selechain,'getcacrd');
    if($nca<=5) {print prcf "$pdbid -9991 -- CACRD data error (NCA $nca) \n";
		 print "$pdbid -9991 -- CACRD data error (NCA $nca) \n\n";next;}
#
  if (! -e "$structuredatadir/$ligandf"){
    $nolig = "true";
    print "Ligand file not found, skipping protein-ligand analysis\n";
  }
  my $num_plg_clu,$plg_bindingsites,$plg_lgcrd;
  if ($nolig eq "true"){
    print "Skipped protein/ligand binding analysis\n";
  } else {
#-- Protein/Ligand binding in X-ray structure
    ($num_plg_clu,$plg_bindingsites,$plg_lgcrd)=
	&get_pdbligand_bindsites($pdbid,$structuredatadir,$workdir,
				 $cacrd,$ligandf,$lgbindcutoff,
				 $lgatomnocutoff,$lgcntatomnocutoff);
    my $calc_info=$num_plg_clu; my $infocalc=$plg_bindingsites;
    if($calc_info < 0) {print prcf "$pdbid $calc_info -- $infocalc\n";
			print "$pdbid $calc_info -- $infocalc\n\n";next;}

    foreach my $plgcid (sort {$a<=>$b} keys %{$plg_bindingsites}){
	print rcdf3 "PLG_BINDINGSITES> $pdbid --- $plgcid -- @{$$plg_bindingsites{$plgcid}}\n";
#	foreach my $lgcrdid (sort {$a<=>$b} keys %{$$plg_lgcrd{$plgcid}}){
#	    print "PLG_LGCTRS> $pdbid--$plgcid--$lgcrdid--@{$$plg_lgcrd{$plgcid}{$lgcrdid}} \n";
    }
#    print "\n"; exit();
#
  } 
#-- DPA binding sites prediction
my $dpaf=$pdbid.'.dpa';
    if(!-e "$dpadatadir/$dpaf"){
	print prcf "$pdbid -2000 -- DPA file do NOT exist\n";
	print "$pdbid -2000 -- DPA file do NOT exist\n";
	next;}
    my ($ndpaclu,$rank_dpaclu,$dpabindingsites,$dpacenters,
	$extremefita,$extremefitb,$serr,$rcorr,$cuttoppercent)=
	    &get_DPA_bindingsites($pdbid,$dpadatadir,$workdir,$cacrd,$dpaf,
				  $toppercent,$dpabindingcutoff);
    my $calc_info=$ndpaclu; my $infocalc=$rank_dpaclu;
    if($calc_info < 0) {print prcf "$pdbid $calc_info -- $infocalc\n";
			print "$pdbid $calc_info -- $infocalc\n\n";next;}

#    print "DPA_EVDFITTING> $pdbid --- $ndpaclu -- $nca -- $extremefita,$extremefitb,$serr,$rcorr \n";
#############
    foreach my $cid (sort {$a<=>$b} keys %$dpabindingsites){
	my $numdpabdsites=scalar(@{$$dpabindingsites{$cid}});
	print rcdf3 "DPA_BINDINGSITES> $pdbid -- $cid -- $numdpabdsites -- @{$$dpabindingsites{$cid}}\n";
#	foreach my $crdid (sort {$a<=>$b} keys %{$$dpacenters{$cid}}){
#	    print rcdf3 "BINDING_DPACRD> $pdbid -- $cid -- $crdid: @{$$dpacenters{$cid}{$crdid}}\n";
#	}
    }    
#    exit();
#############

#
    my $nlgclu,$ndpaclu,$npred,$cmpbindingsites,
	$nullpp,$nullpp2,$ncompp,$prank;
  if ($nolig eq "true"){
    print "Skipping protein-ligand comparison\n";
  } else {
#-- Comparize DPA with X-Ray P/L interactions

    ($nlgclu,$ndpaclu,$npred,$cmp_dpa_plg,$rcd_ncom,$rcd_rank )=
	&comp_xray_dpa($pdbid,$plg_bindingsites,$dpabindingsites,$rank_dpaclu);
    print "$pdbid -1 -- NO DPA_PLG prediction \n" if $npred <=0;
    print "$pdbid $npred -- DPA_PLG prediction -- TOPPERCENT $cuttoppercent\n" if $npred >0;
    print prcf "$pdbid -1 -- NO DPA_PLG prediction \n" if $npred <=0;
    print prcf "$pdbid $npred -- DPA_PLG prediction-- TOPPERCENT $cuttoppercent \n" if $npred >0;

    foreach my $predid (sort {$a<=>$b} keys %{$cmp_dpa_plg}) {
	print "COMP_DPA_PLG> $pdbid -- $predid -- @{$$cmp_dpa_plg{$predid}}\n";
	print rcdf2 "COMP_DPA_PLG> $pdbid -- $predid -- @{$$cmp_dpa_plg{$predid}}\n";
    }
  } 
    printf "%4s %3d %3d %3d TOPP %5.2f  NCOM -- @$rcd_ncom  RANK -- @$rcd_rank\n\n",
    $pdbid,$nlgclu,$ndpaclu,$npred,$cuttoppercent;
    printf rcdf "%4s %3d %3d %3d  %4s %5.2f  %7s %7.3f %7.3f %7.3f %7.3f  NCOM -- @$rcd_ncom  RANK -- @$rcd_rank\n",
    $pdbid,$nlgclu,$ndpaclu,$npred,"TOPP",$cuttoppercent,"FITTING",$extremefita,$extremefitb,$serr,$rcorr;
    print rcdf3 "\n";

    $Nlgclutt+=$nlgclu;$Ndpaclutt+=$ndpaclu;$Npredtt+=$npred;
    $predtt++ if ($npred > 0);
}

my $Rtt=$Ptt=0;
$Rtt=$Npredtt/$Nlgclutt  if($Nlgclutt > 0);
$Ptt=$Npredtt/$Ndpaclutt if($Ndpaclutt > 0);
my $ratio=-1; my $nprot=scalar(@idlist);
$ratio=$predtt/$nprot if($nprot > 0);

print  "\nR-P: $Rtt --- $Ptt\n";
print "Pred: $nprot -- $predtt -- $ratio\n";
print rcdf "\nR-P: $Rtt --- $Ptt\n";
print rcdf "Pred: $nprot -- $predtt -- $ratio\n";


close(rcdf);close(rcdf2);close(rcdf3);close(prcf);
exit();

#
#-----
sub comp_xray_dpa{
    my ($jobid,$ligandbindingsites,$dpabindingsites,$rank_dpaclu)=@_;
    my %cmp_dpa_plg=(); undef %cmp_dpa_plg;

    my $nlgclu=$ndpaclu=0;
    foreach my $cid (sort {$a<=>$b} keys %{$ligandbindingsites}){
	    $nlgclu++;
#	    print "COMP_P_L> $jobid -- $nlgclu -- $cid -- @{$$ligandbindingsites{$cid}}\n";
	}
    foreach my $cid (sort {$a<=>$b} keys %{$dpabindingsites}){
	$ndpaclu++;
#	print "COMP_P_D> $jobid -- $ndpaclu -- $cid -- @{$$dpabindingsites{$cid}}\n";
    }    
###    exit();

#    my @nullpp=();undef @nullpp;my @nullpp2=();undef @nullpp2;
    my @rcd_rank=();undef @rcd_rank;
    my @rcd_ncom=();undef @rcd_ncom;  
    my $ncom=$recall=$precise=$ctdpa_cid='undefine';
    $ncom=$recall=$precise=$ctdpa_cid=-9 if($ndpaclu <=0);

    if($ndpaclu <=0){
	for($i=1;$i<=$nlgclu;$i++){push(@rcd_ncom,-1);push(@rcd_rank,-1)};
	push(@{$cmp_dpa_plg{-1}},$ctdpa_cid,$ncom,$recall,$precise);
	return $nlgclu,$ndpaclu,$npred,\%cmp_dpa_plg,\@rcd_ncom,\@rcd_rank;
    }
    my $npred=0;
    foreach my $plg_cid (sort {$a<=>$b} keys %{$ligandbindingsites}){
	my @xlgbdsites=(); undef @xlgbdsites;
	push(@xlgbdsites,@{$$ligandbindingsites{$plg_cid}});
	($ncom,$recall,$precise,$ctdpa_cid)=
	    &cmp_dpa_single_plgbd(\@xlgbdsites,$dpabindingsites);
	$npred++ if $ncom > 0;	
#	print "COMP_D_L0> $jobid -- $plg_cid -- $ctdpa_cid -- $ncom,$recall,$precise\n";

#	$nullsp=-1;$nullsp2=-1;
#	if($ncom > 0) {    #calc. NULL model
#	    $nlg=scalar(@xlgbdsites);
#	    $ndpa=scalar(@{$$dpabindingsites{$ctdpa_cid}});
#	    $nullsp=&nullmodel($allbindingsites,$ndpa,$nlg,$ncom);
#	    $nullsp2=$nlg/$allbindingsites;
#	}
#	push(@nullpp2,$nullsp2);push(@nullpp,$nullsp);
#
	push(@rcd_ncom,$ncom);
	push(@rcd_rank,@{$$rank_dpaclu{$ctdpa_cid}}[-1]);
	push(@{$cmp_dpa_plg{$npred}},$plg_cid,$ctdpa_cid,
	     @{$$rank_dpaclu{$ctdpa_cid}}[-1],' COMP ',$ncom,$recall,$precise);   #,$nullsp);
	#print "COMP_D_L1> $jobid -- $plg_cid -- $ctdpa_cid -- @{$$rank_dpaclu{$ctdpa_cid}}[-1] ' COMP ' $ncom $recall $precise\n";  #,$nullsp\n";
	#print "COMP_D_L1> $jobid -- $npred -- @{$cmp_dpa_plg{$npred}} \n";
    }
#    print "CMPXD> NULL -- @nullpp  NULL2 -- @nullpp2\n";
#    print "CMPXD> NCOM -- @rcd_ncom\n";
#    foreach my $predid (sort {$a<=>$b} keys %{cmp_dpa_plg}) {
#	print "COMP_DPA_PLG1> $cid -- $predid --  @{$cmp_dpa_plg{$predid}} \n"; }# exit();

    return $nlgclu,$ndpaclu,$npred,\%cmp_dpa_plg,\@rcd_ncom,\@rcd_rank  #,\@nullpp,\@nullpp2;
}
sub cmp_dpa_single_plgbd{    #slgbd: single P/L binding cluster (sites)
    my ($xlgbd0,$sddpabd)=@_;
    my @xlgbd=(); undef @xlgbd; @xlgbd=@{$xlgbd0};
    my $ncom=0;my $recall=0;my $precise=0;my $ctclusterid=-1;
    foreach my $cid (keys %$sddpabd){
	my @dpabd=@{$$sddpabd{$cid}}; #print "@xlgbd \n";	print "@dpabd \n";
	my ($ssncom,$ssrecallLG,$sspreciseDPA)= &cmp_twosets(\@xlgbd,\@dpabd);
#	if($ncom < $ssncom and $precise < $sspreciseDPA){
	if($ncom < $ssncom){
	    $ncom=$ssncom;$ctclusterid=$cid;$recall=$ssrecallLG;$precise=$sspreciseDPA;
	}
    }
    if ($ncom<=0){
	$recall=0;$precise=0;$ctclusterid=-1;
	}
    return $ncom,$recall,$precise,$ctclusterid
}
sub cmp_twosets{
    my ($set0,$set1)=@_;
    my $ncom=0;
    foreach my $il (@$set0){
	foreach my $id (@$set1){
	    if($id == $il) {   
		$ncom++;
		last;
	    }
	}
    }
    my $r0=$r1='undefine';
    my $n0=scalar(@{$set0});my $n1= scalar(@$set1);
    $r0=$ncom/$n0 if $n0 > 0; #recall or precision
    $r1=$ncom/$n1 if $n1 > 0;
    $r0=-1 if($n0 <= 0);
    $r1=-1 if($n1 <= 0); 
    return $ncom,$r0,$r1;
}
sub nullmodel{
    my ($n,$nd,$nl,$nc)=@_;                #$n--total; $nl--ligand; $nd--dpa; $nc--common
    my $n1=$nlg; $n1=$nd if($nd < $nlg);
    my $n2=$n-$nl;
    my $nullpp=0;
    for ($nt=$nc;$nt<=$n1;$nt++){    #n,nd,nl,n1,n2,nt
	$nsp=1;
	for ($i=1;$i<=$nt;$i++){ $nsp=$nsp*($nl-$i+1)/($n-$i+1)};
	for ($i=1;$i<=$nt;$i++){$nsp=$nsp*($nd-$i+1)/($nt-$i+1)};
	if ($nd > $nt) {for ($i=1;$i<=($nd-$nt);$i++){$nsp=$nsp*($n2-$i+1)/($n-$nt-$i+1)} }
	$nullpp+=$nsp;
    }
    return $nullpp;
}

sub bindingstatistics{
    my $ligandbindingsites=@_[0];
    my $dpabindingsites=@_[1];
    my %bindingsitesta=(); undef %bindingsitesta;
    foreach my $df (keys %$ligandbindingsites){
######
#	print "\nBSTA:DOMAIN -- $df \n";
######
	my @lgsitelist=();
	my @dpasitelist=();
	foreach my $role (keys %{$$ligandbindingsites{$df}}) {
	    push (@lgsitelist,$role);
	}
	foreach my $role (keys %{$$dpabindingsites{$df}}) {
	    push (@dpasitelist,$role);
	}
######
#	print "BESTA: LGSTS - @lgsitelist\n";
#	print "BESTA: DPASTS - @dpasitelist \n";
######
	my $ncom=0; my $ratio1=-1; my $ratio2=-1;
	$ncom=cmpsets(\@lgsitelist,\@dpasitelist);
	$ratio1= $ncom/($#lgsitelist+1) if ($#lgsitelist >=0);
	$ratio1= -1.0 if ($#lgsitelist <=-1);
	$ratio2= $ncom/($#dpasitelist+1) if ($#dpasitelist >=0);
	$ratio2= -1.0 if($#dpasitelist<=-1);$ratio2=1 if $ratio2>1;
#	print "BESTA: $ncom $ratio1 $ratio2\n";
	push(@{$bindingsitesta{$df}},$ncom,scalar(@lgsitelist),scalar(@dpasitelist),$ratio1,$ratio2);
    }
    return \%bindingsitesta;
}
sub cmpsets{
    my $lgsitelist=@_[0];
    my $dpasitelist=@_[1];
    my $ncom=0;
    foreach my $il (@$lgsitelist){
	foreach my $id (@$dpasitelist){
	    if($id == $il) {   
#	    if(abs($id-$il) < 2){
		$ncom++;
		last;
	    }
	}
    }
    return $ncom;
}

sub get_DPA_bindingsites{
    my ($jobid,$dpadatadir,$workdir,$cacrd,$dpaf,$toppercent,$dpabindingcutoff)=@_;

#Reading Surface-DPA data
    #system ("cp $dpadatadir/$dpaf $workdir/$dpaf");
    my %surfdpa=(); undef %surfdpa;
    open (sdpaf,"<$dpadatadir/$dpaf");
    my $x=$y=$z='null';
    while (<sdpaf>){  
        #FORTRAN: write(20,'(1x,3F9.3,2x,F12.6,I10)') (hetcrd(k),k=1,3),pert,i
	$x=substr($_,0,10);
	$y=substr($_,10,9);
	$z=substr($_,19,9);
	$pert=substr($_,30,12);
	push(@{$surfdpa{$pert}},$x,$y,$z);
	#print "READDPA> $x\t$y\t$z\t$pert\n";
    }
    close(sdpaf);

    # Write all surface DPA points to PDB if requested
    if ($flag_walldpa == 1) {
        my $alldpapdbf = "$jobid\_alldpa.pdb";
        open(my $alldpa, "> $workdir/$alldpapdbf") || die "Cannot open $workdir/$alldpapdbf\n";
        my $iatom = 0;
        foreach my $pert (sort {$b<=>$a} keys %surfdpa) {
            for (my $k=0; $k<=$#{$surfdpa{$pert}}; $k+=3) {
                $iatom++;
                my $x = $surfdpa{$pert}[$k];
                my $y = $surfdpa{$pert}[$k+1];
                my $z = $surfdpa{$pert}[$k+2];
                printf $alldpa "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s\n",
                    'HETATM', $iatom, ' CA ', ' ', 'DPA', 'A', $iatom, ' ',
                    $x, $y, $z, 1.0, $pert, '', 'SURF';
            }
        }
        close($alldpa);
        print "Wrote all DPA surface points to $workdir/$alldpapdbf\n";
    }
#
#--  Find TOPDPA points and its contact C_ALPHA residue index
    my $cuttoppercent=$toppercent;
increaseperct:
    my ($ntop,$topdpacrd,$extremefita,$extremefitb,$serr,$rcorr)=
	&select_topdpa($jobid,$workdir,\%surfdpa,$ndpasection,$cuttoppercent,$flag_showfit);
    my $infocalc=$topdpacrd if $ntop < 0;
    print "SELE_TOPDPA_ERR> $jobid -- $infocalc\n" if $ntop < 0;
    if ($ntop < 0 && $cuttoppercent > 0.85) {
	$cuttoppercent=$cuttoppercent-0.01;goto increaseperct;
    }
    return -2000, " SELE_TOPDPA_ERR" if $ntop < 0;

#    print "SELE_TOPDPA> $jobid -- $ntop -- $cuttoppercent -- $extremefita,$extremefitb,$serr,$rcorr \n";
#    foreach my $id (sort {$a<=>$b} keys %$topdpacrd){
#	print "SELE_TOPDPA> $jobid -- $id -- @{$$topdpacrd{$id}}\n";
#    } #exit();

# 
#-- Clustering TOP_DPA points
    my $topdpaclusterpdbf='';
    $topdpaclusterpdbf=$jobid.'_topdpacluster.pdb' if ($flag_wclusterpdb == 1);
    my $topdpacluster=OPTICS($topdpacrd,$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa,$topdpaclusterpdbf);

###    print "SELE_TOPDPA> $ntop,$epislon_dpa,$MinPts_dpa,$epislon_cluster_dpa\n";
#    foreach $cid (sort {$a<=$b} keys %$topdpacluster){
###	print "TOPDPA_CLU> $jobid -- $cid --  @{$$topdpacluster{$cid}}\n";}    exit();

#
#-- Protein/DPA_cluster binding sites and Ranking
    my ($dpabindingsites,$dpacenters)=&getbindingdata_with_ClusterPOINTS($topdpacluster,$topdpacrd,$cacrd,$dpabindingcutoff);
    $ndpaclu=0;
    foreach my $cluster_id (sort {$a<=>$b} keys %$dpabindingsites){ 
	$ndpaclu++ if $cluster_id>0;} 
    if ( $ndpaclu <=0 && $cuttoppercent > 0.85) {  
	if ($adjcutperct eq "true"){
	    $cuttoppercent=$cuttoppercent-0.01;
	    print prcf "$jobid -- ADJUST TOP_DPA_PERCENTCUT: $cuttoppercent\n";
	    print "$jobid -- ADJUST TOP_DPA_PERCENTCUT: $cuttoppercent\n";
	    goto increaseperct;}}
    return -1500," NO eff_TOP_DPA_CLUSTER_ERR" if $ndpaclu <= 0;

    my $rank_dpaclu=&rank_of_cluster($ndpaclu,$topdpacluster,$topdpacrd);
#########vvv  
#    foreach my $cid (sort {$a<=>$b} keys %$topdpacluster) {
#	print "TOPDPA_CLU> $cid -- @{$$topdpacluster{$cid}}\n";}
#    foreach my $cid (sort {$a<=>$b} keys %$dpabindingsites){
#	print "DPACLU_BINDINGSITES> $jobid -- $ndpaclu -- $cid -- @{$$dpabindingsites{$cid}}\n";
#	foreach my $crdid (sort {$a<=>$b} keys %{$$dpacenters{$cid}}){
#	    print "DPACLU_DPACRD> $jobid -- $cid -- $crdid: @{$$dpacenters{$cid}{$crdid}}\n";
#	} 
#    } 
#    foreach $cid (sort {$a<=>$b} keys %$rank_dpaclu){
#	print "DPACLU_RANK> $jobid -- $cid -- @{$$rank_dpaclu{$cid}}\n" if $cid > 0;}
#    exit();
#########^^^
#
    return $ndpaclu,$rank_dpaclu,$dpabindingsites,$dpacenters,$extremefita,$extremefitb,$serr,$rcorr,$cuttoppercent; 
}

sub select_topdpa{
    my ($jobid,$workdir,$sdpa,$nsec,$cutperct,$flag_showfit)=@_;
#
#-- Extreme Value Fitting
    my $sdpalistf=$jobid.'_dpa.list';
    my $fitted_file=$jobid.'_EVD.fit';
    my $extremefitinputf=$jobid.'_EVD.fit_inp';
    my $efitout=$jobid.'_EVD.fit_out';

    open(sdpalist,"> $workdir/$sdpalistf")||die "Cannot open $workdir/$sdpalistf";
    foreach my $pert (keys %$sdpa){
	for ($k=0;$k<=$#{$$sdpa{$pert}};$k+=3){print sdpalist "  $pert\n";}
    }close(sdpalist);
    open(extremefitinput,"> $workdir/$extremefitinputf");
    print extremefitinput " \"$workdir/$sdpalistf\" $nsec $cutperct $flag_showfit \n";
    print extremefitinput " \"$workdir/$fitted_file\" \n";
    close(extremefitinput);
    system ("$exedir/extremefit < $workdir/$extremefitinputf  > $workdir/$efitout");
#
#-- Determine DPA threshold
    my $dpa_threshold='undefine';my $serr=9999.;my $rcorr=0.;
    open(tmpf,"< $workdir/$efitout"); my $line=0;my $extremefita = my $extremefitb ='null';
    while(<tmpf>){ 
	my @tmp=split /\s+/,$_;$line++;#print 'READFIT> ',$_;
	if($line==1 && @tmp[0] eq ''){$dpa_threshold=@tmp[1];$serr=@tmp[2];$rcorr=@tmp[3];};
	if($line==1 && @tmp[0] ne ''){$dpa_threshold=@tmp[0];$serr=@tmp[1];$rcorr=@tmp[2];};
	if($line==2 && @tmp[0] eq ''){ $extremefita=@tmp[1];$extremefitb=@tmp[2];last;};
	if($line==2 && @tmp[0] ne ''){ $extremefita=@tmp[0];$extremefitb=@tmp[1];last;};
    }
    close(tmpf);
    print "$jobid calc. segment damp in EVD_fitting\n" if($dpa_threshold eq 'undefine');
    print prcf "$jobid calc. segment damp in EVD_fitting\n" if($dpa_threshold eq 'undefine');
    system ("rm -f $workdir/$extremefitinputf $workdir/$efitout $workdir/$sdpalistf"); 
    #print "FITTING> $dpa_threshold -- $extremefita -- $extremefitb -- $serr -- $rcorr\n";
    return -1, "sele_topDPA error or EVD fitting error" if($dpa_threshold eq 'undefine');
#
#-- Slect topDPA points
    my $ntop=-1; my %topdpa=(); undef %topdpa;
    foreach my $pert (sort {$b<=>$a} keys %$sdpa){
	last if $pert < $dpa_threshold;
	for ($k=0;$k<=$#{$$sdpa{$pert}};$k+=3){
	    $ntop++;
	    push(@{${topdpa{$ntop}}},@{$$sdpa{$pert}}[k],@{$$sdpa{$pert}}[k+1],
		 @{$$sdpa{$pert}}[k+2],$ntop,'TOP',$pert);
	}}
    return -1,"Nothing selected" if $ntop < 0;
###    foreach my $id (sort {$a<=>$b} keys %topdpa){
#	print "SELE_TOPDP> $id -- @{$topdpa{$id}}  \n";
###    } print "\n";
    return $ntop,\%topdpa,$extremefita,$extremefitb,$serr,$rcorr;
}

################################################################
sub get_pdbligand_bindsites{
    my ($jobid,$datadir,$workdir,$cacrd,$ligandf,$lgbindcutoff,$lgatomnocutoff,$lgcntatomnocutoff)=@_;
#
#-- Ligand CRD
    return -8999," Not detect LIGAND file" if(! -e "$datadir/$ligandf");
    my ($nlg,$lgcrd)=&getligandcrd($datadir,$workdir,$ligandf,'getNOHcrd');
    return -8999," Ligand atoms are too small ($nlg)" if $nlg<$lgatomnocutoff;
#    foreach $hid (sort {$a<=>$b} keys %$lgcrd){
#    print "LIGAND_CRD> $hid -- @{$$lgcrd{$hid}}\n";}exit();
#    &writepdbf($workdir,'lgatm.pdb',$lgcrd,'9');
# 
#-- Clustering LIGAND atoms
    my $ligandclusterpdbf='';
    $ligandclusterpdbf=$jobid.'_ligandcluster.pdb' if ($flag_wclusterpdb == 1);
    my $lgcluster=OPTICS($lgcrd,5,2,5,$ligandclusterpdbf);

#    foreach (sort {$a<=$b} keys %$lgcluster){
#	print "LIGAND_CLU> $jobid -- $_ --  @{$$lgcluster{$_}}\n";}    exit();

#
#-- Protein/Ligand binding sites
    my ($plgcts,$plgctrs)=&getbindingdata_with_ClusterPOINTS($lgcluster,$lgcrd,$cacrd,$lgbindcutoff);		

    my %ligandbindingsites=(); undef %ligandbindingsites; #output  HASH: { cluster_id -- array of resid}
    my %bindingligandcrd=();   undef %bindingligandcrd;   #output  HASH of HASH: { cluster_id -- ligandATOM_ID -- crd}  
    $binding_cluster=-1;
    foreach my $cid (sort {$a<=>$b} keys %$plgcts){
	next if($#{$$plgcts{$cid}} < $lgcntatomnocutoff);   #Ignore small lgcnt portions
	push (@{$ligandbindingsites{++$binding_cluster}},@{$$plgcts{$cid}});
	foreach my $rol (sort {$a<=>$b} keys %{$$plgctrs{$cid}}){
	    push (@{$bindingligandcrd{$binding_cluster}{$rol}},@{$$plgctrs{$cid}{$rol}});
	}
    }
##
######################
#    foreach my $cid (sort {$a<=>$b} keys %$plgcts){
#	print "PLI_CTS> $jobid -- $lgbindcutoff || $cid -- @{$$plgcts{$cid}}\n";
#    } print "\n";
#   foreach my $cid (sort {$a<=>$b} keys %$plgctrs){
#       foreach my $rol (sort {$a<=>$b} keys %{$$plgctrs{$cid}}){
#	    print "PLI_CTRS> $jobid -- $cid -- $rol -- @{$$plgctrs{$cid}{$rol}}\n";
#	}
#    } print "\n";
#    foreach my $lgcid (sort {$a<=>$b} keys %{ligandbindingsites}){
#	print "GLG_CTS_> $pdbid -- $lgcid  @{$ligandbindingsites{$lgcid}}\n";
#	foreach my $hcrdid (sort {$a<=>$b} keys %{$bindingligandcrd{$lgcid}}){
#	    print "GLG_CTRS_> $pdbid -- $lgcid -- $hcrdid -- @{$bindingligandcrd{$lgcid}{$hcrdid}} \n";
#	}
#    } print "\n"; #exit();
#######################
#
    return -5999," NO P/L interactions!" if $binding_cluster <0;
    return $binding_cluster,\%ligandbindingsites,\%bindingligandcrd if $binding_cluster >=0;
}
sub getligandcrd{
    my ($datadir,$outputdir,$ligandf,$flag_readpdb,$recordf)=@_;
    my $nlg=-1; #number of ligand atoms
    my %lgcrd=(); my%lgcrd_detail=(); undef %lgcrd; undef %lgcrd_detail;

    my $atom=$serial=$atomname=$altLoc=$resName=$chainID=$resSeq='';
    my $iCode=$x=$y=$z=$occupancy=$tempFactor=$segID=$element=$charge='';
    my $resSeq0=$iCode0=$serial0='null';
    $pdbformat='%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s';
    open(wcrd,"> $outputdir/$recordf") if($recordf ne '');           #wcrd -- writing crd
    open(lgf,"< $datadir/$ligandf");
    while (<lgf>){
	$atom=substr($_,0,6);
	$serial=substr($_,6,5);
	$atomname=substr($_,12,4);
	$altLoc=substr($_,16,1);
	$resName=substr($_,17,3);
	$chainID=substr($_,21,1);
	$resSeq=substr($_,22,4);
	$iCode=substr($_,26,1);
	$x=substr($_,30,8);
	$y=substr($_,38,8);
	$z=substr($_,46,8);
	$occupancy=substr($_,54,6);
	$tempFactor=substr($_,60,6);
	$segID=substr($_,72,4);
	$element=substr($_,76,2);
	$charge=substr($_,78,2);
	next if($atom ne "ATOM  " and $atom ne 'HETATM');
#
#--- ALL ATOM
	if( $flag_readpdb eq 'getall' or $flag_readpdb eq 'getallcrd'){
	    if($resSeq ne $resSeq0  or $iCode ne $iCode0 or $serial ne $serial0) {
		$resSeq0=$resSeq;$iCode0=$iCode;$serial0=$serial;$nlg++;
		push (@{$lgcrd{$nlg}},$x,$y,$z,$resSeq,$resName);
		print wcrd "\t$x\t$y\t$z\t$resSeq\t$resName \n" if($recordf ne '');}
	}
	elsif($flag_readpdb eq 'getNOH' or $flag_readpdb eq 'getNOHcrd'){   #Not Hydrogen atom
	    if($resSeq ne $resSeq0  or $iCode ne $iCode0 or $serial ne $serial0) {
		$resSeq0=$resSeq;$iCode0=$iCode;$serial0=$serial;
		my $hid=substr($atomname,1,1);
		if($hid ne 'H' and $hid ne 'h' and $element ne ' H' and $element ne ' h' and $element ne 'H ' and $element ne 'h '){
		    $nlg++;
		    push (@{$lgcrd{$nlg}},$x,$y,$z,$resSeq,$resName);
		    print wcrd "\t$x\t$y\t$z\t$resSeq\t$resName \n" if($recordf ne '');}
	    }
	}
    }
    close(wcrd) if($recordf ne '');
    return $nlg,\%lgcrd;
}




















#########################################################
#########################################################
#########################################################
sub timetitle{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime; #localtime(time);
    $year = $year - 100; # Handle century
    $mon++;
    $hour=$hour-6;
    if (length($year) < 2) { $year = "0" . $year;} # Make the years all two digits
    if (length($mon) < 2) { $mon = "0" . $mon;} # Make the months all two digits 
    if (length($mday) < 2) { $mday = "0" . $mday;} # Make the dates all two digits
    my $timetile= "$year-$mon-$mday $hour:$min";
    return $timetile;
}
sub writepdbf{
    my ($outputdir,$wpdbfilename,$crdhash,$chainid)=@_;
    open(wp,"> $outputdir/$wpdbfilename");
    my $lid=0;
    foreach $hid (sort {$a<=>$b} keys %$crdhash) {
	my $x=@{$$crdhash{$hid}}[0];
	my $y=@{$$crdhash{$hid}}[1];
	my $z=@{$$crdhash{$hid}}[2];
	my $resSeq=@{$$crdhash{$hid}}[3];
	my $resName=@{$$crdhash{$hid}}[4];
	$lid++;
	printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n",
	'HETATM',$lid,' CA ',' ',$resName,$chainid,$resSeq,' ',$x,$y,$z,1.0,10.0,'','MING','','';
#	printf  "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 'ATOM  ',
#	$lid,' CA ',' ',$resName,' ',$resSeq,' ',$x,$y,$z,1.0,10.0,'','MING','','';
    }
    close(wp);
}
sub writepdbf2{
    my $wpdbfilename=@_[0];
    my $crdhash=@_[1];
    my $chainid=@_[2];
    my $thred = @_[3];
    open(wp,"> $wpdbfilename");
    $lid=0;$resName='ALA';
    foreach $hid (sort {$b<=>$a} keys %$crdhash) {
	my $x=@{$$crdhash{$hid}}[0];
	my $y=@{$$crdhash{$hid}}[1];
	my $z=@{$$crdhash{$hid}}[2];
	my $pert=@{$$crdhash{$hid}}[3];
	next if $pert < $thred;
#	my $resName=@{$$crdhash{$hid}}[4];
	$lid++;$resSeq=$lid;
	printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 'HETATM',
	$lid,' CA ',' ',$resName,$chainid,$resSeq,' ',$x,$y,$z,1.0,$pert,'','MING','','';
#	printf  "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 'ATOM  ',
#	$lid,' CA ',' ',$resName,' ',$resSeq,' ',$x,$y,$z,1.0,10.0,'','MING','','';
    }
    close(wp);
}





#########################################
############# OPTICS ####################
#########################################

############### CLUSTER OF TOPDPA POINTS #######################
sub OPTICS{
    ($crddata,$epislon,$MinPts,$epislon_cluster,$flag_wpdb)=@_;               #Global PARAMETERS in Cluster
    $nbMinPts=$MinPts-1;
    $epislons=$epislon**2;$epislon_clusters=$epislon_cluster**2;
    $for_core_distanceSq='core_distanceSq';
    $if_processed='if_processed';
    $for_reach_distSq='reachable_distanceSq';
    %oset=();%infoset=();%seedlist=();$nlist=0;%listseed=();@orderlist=();      #Global LIST of database 
    undef %oset;undef %infoset;undef %seedlist; undef %listseed;undef @orderlist;
    my $nop=-1; foreach $i ( keys %$crddata){push (@{$oset{$i}},@{$$crddata{$i}}); $nop++;}
#########
#    foreach (sort {$a<=>$b} keys %oset){ print "CHECK_OSET>$_:    @{$oset{$_}}\n"; }
#########

    for(my $i=0;$i<=$nop;$i++){
	next if ($infoset{$i}{$if_processed} eq 'yes');
	&ExpandClusterOrder($i);
    }
    my $optics_cluster=OPTICS_CLUSTER();
    &pdb_showCluster($flag_wpdb,$optics_cluster) if($flag_wpdb ne '');

#########
#    foreach (keys %$optics_cluster){	print "OPTICS> $_  @{$$optics_cluster{$_}}\n";}
#    print "ORDERLIST> @orderlist\n"; 
#########

    return $optics_cluster;
}
sub OPTICS_CLUSTER{
    my %cluster=();undef %cluster;my $clusterid=0; #0 for NOISE, 1,2,3... for nontrivial cluters
    
    foreach my $i (@orderlist){
	$rcsq=$infoset{$i}{$for_core_distanceSq};
	$rdsq=$infoset{$i}{$for_reach_distSq} ;
	if($rdsq > $epislon_clusters or $rdsq eq 'undef'){
	    if($rcsq<=$epislon_clusters and $rcsq ne 'undef' ){
		$clusterid++;
		push(@{$cluster{$clusterid}},$i);}
	    else{
		push(@{$cluster{0}},$i);}
	}
	else{
	    push(@{$cluster{$clusterid}},$i);}
    }
    return \%cluster;
}
sub ExpandClusterOrder{
    my $id=@_[0];
    my ($nb,$neighbor)=Neighbor($id);
    $infoset{$id}{$if_processed}='yes';
    $infoset{$id}{$for_core_distanceSq}=CoreDistance($neighbor);
    push (@orderlist,$id); #print "ADDING ($id) to orderlist\n";
    $infoset{$id}{$for_reach_distSq}='undef';
    return  if($infoset{$id}{$for_core_distanceSq} eq 'undef');   #NOT A core-object
    &OrderSeedUpdata($neighbor,$id);
    return if($nlist <=0);
    while ($nlist >=1) {
	($first_rds,$first_listid)=NextSeedlist(); 
	splice(@{$seedlist{$first_rds}},0,1);$nlist--;
	delete($listseed{$first_listid});
	if(scalar(@{$seedlist{$first_rds}}) eq 0){
	    delete($seedlist{$first_rds});}
	($nb,$neighbor)=Neighbor($first_listid);
	$infoset{$first_listid}{$if_processed}='yes';
	$infoset{$first_listid}{$for_core_distanceSq}=CoreDistance($neighbor);
	push (@orderlist,$first_listid); #print "whileloop)ADDING ($first_listid) with CORE $infoset{$first_listid}{$for_core_distanceSq} to orderlist\n" if($id ne '');

	next if($infoset{$first_listid}{$for_core_distanceSq} eq 'undef');
	&OrderSeedUpdata($neighbor,$first_listid);
    }
}
sub CoreDistance{
    my $neighbor=@_[0];
    my $nnb=0; my $rcores='undef';
    foreach my $rs (sort {$a<=>$b} keys %$neighbor){
	$nnb+=scalar(@{$$neighbor{$rs}});
	if($nnb >= $nbMinPts){      #counting the object itself
	    $rcores=$rs;
	    last};
    }
    return $rcores;
}
sub Neighbor{     #need: %oset,$epislons
    my $id=@_[0];
    my %neighbor=(); undef %neighbor;
    my $x0,$y0,$z0,$x,$y,$z,$dx,$dy,$dz;my $nb=0;
    $x0=@{$oset{$id}}[0];$y0=@{$oset{$id}}[1];$z0=@{$oset{$id}}[2];
    foreach my $i (keys %oset){
	next if($i eq $id);
	$x=@{$oset{$i}}[0];$y=@{$oset{$i}}[1];$z=@{$oset{$i}}[2];
	$dx=$x0-$x;$dy=$y0-$y;$dz=$z0-$z;
	next if($dx>$epislon or $dy>$epislon or $dz>$epislon or $dx<-$epislon or $dy<-$epislon or $dz<-$epislon);
	$ds=$dx**2+$dy**2+$dz**2;
	next if($ds > $epislons);
	push (@{$neighbor{$ds}},$i);
	$nb++;
    }
    return $nb,\%neighbor;
}
sub OrderSeedUpdata{
    my $neighbor=@_[0];my $id=@_[1];
    my $reach_distSq,$old_rds;
    my $core_distSq=$infoset{$id}{$for_core_distanceSq};  #$core_distSq must PREDEFINED

    foreach my $ds (sort {$a<=>$b} keys %$neighbor){
	foreach my $nbid (@{$$neighbor{$ds}}){
	    $reach_distSq=$core_distSq;$reach_distSq=$ds if($ds > $core_distSq);
	    if ($infoset{$nbid}{$if_processed} eq 'yes'){
		next if ($infoset{$nbid}{$for_reach_distSq} ne 'undef' or $infoset{$nbid}{$for_core_distanceSq} ne 'undef');
		for(my $tid=0;$tid<=$#orderlist;$tid++){
		    if(@orderlist[$tid] eq $nbid){
			splice(@orderlist,$tid,1);last;}}
	    }
	    if($listseed{$nbid} eq ''){
		$listseed{$nbid}=$reach_distSq;
		$infoset{$nbid}{$for_reach_distSq}=$reach_distSq;    #renew reachable_distance
		push(@{$seedlist{$reach_distSq}},$nbid);
		$nlist++;
	    }
	    elsif($reach_distSq <$listseed{$nbid}){
		$old_rds=$listseed{$nbid};
		$listseed{$nbid}=$reach_distSq;
		$infoset{$nbid}{$for_reach_distSq}=$reach_distSq;
		my $tn=-1;
		for(my $tid=0;$tid<=$#{$seedlist{$old_rds}};$tid++){
		    if(@{$seedlist{$old_rds}}[$tid] eq $nbid){
			splice(@{$seedlist{$old_rds}},$tid,1);$nlist--;last;}}
		if(scalar(@{$seedlist{$old_rds}}) eq 0){
		    delete($seedlist{$old_rds});}
		push(@{$seedlist{$reach_distSq}},$nbid);$nlist++; #adding $nbid with new Rds
	    }
	}
    }
    return;
}
sub NextSeedlist{
    my @rds=(sort {$a<=>$b} keys %seedlist);
    my $first_rds=@rds[0];
    my $first_listid=@{$seedlist{@rds[0]}}[0];
    return $first_rds,$first_listid;
}
sub pdb_showCluster{          
    my $wpdbfilename=@_[0];
    my $cluster=@_[1];   
    my @chainid_list=(O..Z);my @resName_list=('CLU','DLU','XLU','ZLU');
    open(wp,"> $workdir/$wpdbfilename");
    my $iatom=0;
    foreach $cid (sort {$a<=>$b} keys %$cluster) {
	$chainid=X if($cid eq 0);
	$chainid=@chainid_list[($cid-1)%scalar(@chainid_list)] if($cid >=1);
	$fold_id=int(($cid-1)/scalar(@chainid_list));
	$resname=@resName_list[$fold_id];
	$resname=@resName_list[$fold_id%scalar(@resName_list)] if($fold_id >= scalar(@resName_list));;
	$resname='NOS' if($cid eq 0);
	my $ia_chain=0;
	foreach (my $i=0;$i<=$#{$$cluster{$cid}};$i++){
	    my $id=@{$$cluster{$cid}}[$i];$pid=$id;
	    my $x=@{$oset{$id}}[0]; my $y=@{$oset{$id}}[1]; my $z=@{$oset{$id}}[2];
	    my $resSeq=@{$oset{$id}}[3];
	    my $resName=@{$oset{$id}}[4];
	    my $value=10.00;
	    my $value=@{$oset{$id}}[5] if @{$oset{$id}}[5] ne '';
	    $iatom++;
	    $ia_chain++;

#	    printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n", 'HETATM',
#	    $iatom,' CA ',' ',$resname,$chainid,$ia_chain,' ',$x,$y,$z,1.0,10.0,'','MING','_C','LA';

	    printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%4d\n", 'HETATM',
	    $iatom,' CA ',' ',$resName,$chainid,$resSeq,' ',$x,$y,$z,1.0,$value,'','MING',$pid;

#	    printf wp "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%4d\n", 'HETATM',
#	    $iatom,' CA ',' ',$resName,$chainid,$ia_chain,' ',$x,$y,$z,1.0,$value,'','MING',$pid;

	    }
    }
    close(wp);
}
############### END OF CLUSTER PROGRAM #######################
sub read_pdb{
    my ($pdbdatadir,$pdbf,$selechain,$flag_read,$recordf,$outputdir)=@_;
    $selechain='FirstChain' if ($selechain eq '' || $selechain eq 'undefine');

    my $nca=$nhet=$natom=$nenv=0;
    my %cacrd=%hetcrd=%envcrd=@envid=(); undef %cacrd;undef %hetcrd;undef %envcrd; undef @envid;

    my $atom=$serial=$atomname=$altLoc=$resName=$chainID=$chainIDw=$resSeq='';
    my $iCode=$x=$y=$z=$occupancy=$tempFactor=$segID=$element=$charge='';
    my $resSeq0=$iCode0='null'; my $atmtype='null';
    my @revisedRES=('ACE','CGU','MSE','CME','CSS','KCX','TRO','SEP',' FE','OXY');

    open(wcrd,"> $outputdir/$recordf") if($recordf ne '' and $recordf ne 'none');
    open(pf,"<  $pdbdatadir/$pdbf");
    while (<pf>){
	$atom=substr($_,0,6);
	$serial=substr($_,6,5);
	$atomname=substr($_,12,4); $atmtype=substr($_,12,1); $atmtype=substr($_,13,1) if $atmtype=' ';
	$altLoc=substr($_,16,1);
	$resName=substr($_,17,3);
	$chainID=substr($_,21,1);
	$resSeq=substr($_,22,4);
	$iCode=substr($_,26,1);
	$x=substr($_,30,8);
	$y=substr($_,38,8);
	$z=substr($_,46,8);
	$occupancy=substr($_,54,6);
	$tempFactor=substr($_,60,6);
	$segID=substr($_,72,4);
	$element=substr($_,76,2);
	$charge=substr($_,78,2);
	next if($atom ne "ATOM  " and $atom ne 'HETATM');
	$selechain=$chainID if $selechain eq 'FirstChain';
	$chainIDw=$chainID; $chainIDw='0' if $chainIDw eq ' ';

	if(($atom eq "ATOM  " or grep (/$resName/, @revisedRES)) and $resName ne 'ACE'  and  $atomname eq ' CA ' and ($flag_read eq 'getnca' or $flag_read eq 'getcacrd')){
	    if($resSeq ne $resSeq0 or $iCode ne $iCode0) {
		$resSeq0=$resSeq;
		$iCode0=$iCode;
		if ($selechain eq $chainID || $selechain eq 'ALL'|| $selechain eq 'all' || $selechain eq 'All'){
		    $nca++   if($flag_read eq 'getnca' or $flag_read eq 'getcacrd');
		    push (@{$cacrd{$nca}},$x,$y,$z,$resSeq,$tempFactor,$resName) if($flag_read eq 'getcacrd'); 
		    print wcrd "\t$x\t$y\t$z\t$resSeq\t$tempFactor\t$chainIDw\t$resName\n" if($recordf ne '' and $recordf ne 'null');
		}
	    }
	}
	if(($atom eq 'HETATM' and ( $hetatm_name eq '***' or  $resName eq $hetatm_name) and !grep (/$resName/, @revisedRES) and $resName ne 'HOH') and ($flag_read eq 'getnhet' or $flag_read eq 'gethetcrd')){
	    if ($selechain eq $chainID || $selechain eq 'ALL'|| $selechain eq 'all' || $selechain eq 'All'){
		$nhet++   if($flag_read eq 'getnhet' or $flag_read eq 'gethetcrd');
		push (@{$hetcrd{$nhet}},$x,$y,$z,$resSeq,$resName) if($flag_read eq 'gethetcrd' or $flag_read eq 'getall');
		print  wcrd "\t$x\t$y\t$z\t$resSeq\t$tempFactor\t$chainIDw\t$resName\n" if($recordf ne '' and $recordf ne 'none');
	    }
	}
	if(($atom eq "ATOM  " or grep (/$resName/, @revisedRES)) and $atomname ne ' CA ' and $atmtype ne 'H' and ($flag_read eq 'getnheavyenv' or $flag_read eq 'getheavyenvcrd')){
	    if ($selechain eq $chainID || $selechain eq 'ALL'|| $selechain eq 'all' || $selechain eq 'All'){ 
		$nenv++;
		push(@envid,$serial);
		push (@{$envcrd{$nenv}},$x,$y,$z,$resSeq,$resName) if($flag_read eq 'getheavyenvcrd');
		print wcrd "\t$x\t$y\t$z\t$resSeq\t$atom\t$chainIDw\t$atmtype\t$tempFactor\n" if($recordf ne '' and $recordf ne 'null');
	    }
	}
	if(($atom eq "ATOM  " or grep (/$resName/, @revisedRES)) and $atomname ne ' CA 'and ($flag_read eq 'getnenv' or $flag_read eq 'getenvcrd')){
	    if ($selechain eq $chainID || $selechain eq 'ALL'|| $selechain eq 'all' || $selechain eq 'All'){ 
		$nenv++;
		push(@envid,$serial);
		push (@{$envcrd{$nenv}},$x,$y,$z,$resSeq,$resName) if($flag_read eq 'getenvcrd');
		print wcrd "\t$x\t$y\t$z\t$resSeq\t$atom\t$chainIDw\t$atmtype\t$tempFactor\n" if($recordf ne '' and $recordf ne 'null');
	    }
	}
	if($flag_read eq 'getpeptidecrd' and ($atom eq "ATOM  " or grep (/$resName/, @revisedRES))){
	    $natom++;
	    print wcrd "\t$x\t$y\t$z\t$resSeq\t$resName\t$chainIDw\n";
	}
	if($flag_read eq 'getpeptidepdb' and ($atom eq "ATOM  " or grep (/$resName/, @revisedRES))){
	    print wcrd "$_";
	    $natom++;
	}
    }
    close(pf);
    close(wcrd) if($recordf ne '' and $recordf ne 'none');
    return $nca  if($flag_read eq 'getnca');
    return $nca,\%cacrd  if($flag_read eq 'getcacrd');
    return $nhet  if($flag_read eq 'getnhet');
    return $nhet,\%hetcrd  if($flag_read eq 'gethetcrd');
    return $nenv  if($flag_read eq 'getnenv' or $flag_read eq 'getheavynenv');
    return $nenv,\%envcrd  if($flag_read eq 'getenvcrd' or $flag_read eq 'getheavyenvcrd');
    return $natom if($flag_read eq 'getpeptidecrd' or $flag_read eq 'getpeptidepdb');
}
sub getbindingdata_with_ClusterPOINTS{
    my ($targetcluster,$targetcrd,$cacrd,$bcut)=@_;
#    my $targetcluster=@_[0];
#    my $targetcrd=@_[1];
#    my $cacrd=@_[2];
#    my $bcut=@_[3];
    my %targetcts=(); undef %targetcts; my %targetctrs=(); undef %targetctrs;
#    foreach my $cid (sort {$a<=>$b} keys %$targetcluster) {print "GCPCTS_CHECK_CLU> $cid -- @{$$targetcluster{$cid}}\n"; last if $cid >3; }
#    foreach my $cid (sort {$a<=>$b} keys %$targetcrd){print "GCPCTS_CHECK_TCRD>  $cid   @{$$targetcrd{$cid}}\n"; last if $cid >3;}
#    foreach my $cid (sort {$a<=>$b} keys %$cacrd) {print "GCPCTS_CHECK_CA> $cid -- @{$$cacrd{$cid}}\n"; last if $cid >3;}
#    foreach my $cid (sort {$a<=>$b} keys %$targetcluster) {print "GCPCTS_CHECK_CLU> $cid -- @{$$targetcluster{$cid}}\n";}
    foreach my $targetclusterster_id (sort {$a<=>$b} keys %$targetcluster){
	next if($targetclusterster_id eq 0); #print "GTARGETCTS_CHECK> $targetclusterster_id -- @{$$targetcluster{$targetclusterster_id}}\n";
	my %targetcenters=(); undef %targetcenters; my $np=-1;
	foreach $crd_id (@{$$targetcluster{$targetclusterster_id}}){
	    push (@{$targetcenters{++$np}},@{$$targetcrd{$crd_id}}); #print "GTARGETCTS> $targetclusterster_id -- $crd_id -- $np -->@{$$topcrd{$crd_id}}\n";
	    push (@{$targetctrs{$targetclusterster_id}{$crd_id}},@{$$targetcrd{$crd_id}});
	}
	my $targetcontacts=get_Ligand_contacts_Calpha($cacrd,\%targetcenters,$bcut); #print "GTARGETCTS> $targetclusterster_id --  @{$targetcontacts}   <<---\n";
	$targetcts{$targetclusterster_id}=$targetcontacts;
    }
#########vvv
#    foreach $did (sort {$a<=>$b} keys %targetcts){
#	print "GTARGETCTS_> $did  - @{$targetcts{$did}}\n";
#    } 
#    foreach $cid (sort {$a<=>$b} keys %targetctrs){
#	foreach $rol (sort {$a<=>$b} keys %{$targetctrs{$cid}}){
#	    print "GTARGETCTRS_> $cid -- $rol --- @{$targetctrs{$cid}{$rol}}\n";
#	}
#    } 
#########^^^
    return \%targetcts,\%targetctrs;
}
sub get_Ligand_contacts_Calpha{
    my $cacrd=@_[0];                    #Hash 1: %$cacrd   
    my $hetcrd=@_[1];                   #Hash 2: %$hetcrd
    my $bindingcutoff=@_[2];
    my $bindingcutoffsq=$bindingcutoff**2;
#########
#    foreach (sort{$a<=>$b} keys %$hetcrd){
#	print "HETCRD> $bindingcutoff || $_ --  @{$$hetcrd{$_}}\n"
#	}
#    foreach (sort{$a<=>$b} keys %$cacrd){
#	print "CA_CRD> $bindingcutoff || $_ --  @{$$cacrd{$_}}\n"
#	}
#    exit();
#########
    my @contacts=();                    #output -- Recording index in cacrd and the distance
    my @d1=();my $dtmp=''; my $dissq='';

    foreach my $caid (keys %$cacrd){              #For each CA searching ALL Ligand atoms
	foreach $hetid (keys %$hetcrd){
	    @d1=();                              #tmp displacement between Ligand and CA
	    for (my $k=0;$k<=2;$k++){
		$dtmp=@{$$cacrd{$caid}}[$k]-@{$$hetcrd{$hetid}}[$k];
		last if($dtmp > $bindingcutoff or $dtmp < -$bindingcutoff);
		push(@d1,$dtmp);
	    }
	    next if($#d1 <2);                     # at least 1 of 3 components of displacement > cutoff
	    $dissq=(@d1[0]*@d1[0]+@d1[1]*@d1[1]+@d1[2]*@d1[2]);
	    if($dissq<$bindingcutoffsq){          #find the contact CA
		push (@contacts,$caid);   #   print "GCTS> $bindingcutoff -- $caid\n";
		last;                     #stop search when any ligand atom found within the cutoff
	    }
	}
    }
######
#    print "GETCONTACT: $bindingcutoff -- @contacts\n"; 
######
    return  \@contacts;
}

sub rank_of_cluster{
    my ($nclu,$crdclu,$crd)=@_;
    my %clurank=(); undef %clurank;
   
    push (@{$clurank{-1}},-1,-1,-1);
    if ($nclu<=0){
	push(@{$clurank{0}},-1,-1);
	return \%clurank;}
    if ($nclu==1){
	push(@{$clurank{1}},1,1);	
	return \%clurank;}

    my $cluminino=9999;
    foreach my $cid  (keys %$crdclu){
	next if $cid eq 0;
	my $tmpn=scalar(@{$$crdclu{$cid}});
	$cluminino=$tmpn if $tmpn < $cluminino;
    }
    my @pertavgall=();undef @pertavgall;
    foreach my $cid (sort {$a<=>$b} keys %$crdclu){    #DPADATA ordered as "pert" values
	next if $cid <=0;
	my $perttt=0.; my $pertavg=-1; my %tmpcrd=(); undef %tmpcrd;
	my $ncrd=0;                                   #take first $cluminino data
	foreach my $crdid (@{$$crdclu{$cid}}){
	    $perttt+=@{$$crd{$crdid}}[-1]; $ncrd++;last if $ncrd > $cluminino; #take first $cluminino data
	    push (@{$tmpcrd{$crdid}},@{$$crd{$crdid}});
	    #print "RANK> $crdid -- $ncrd -- @{$tmpcrd{$crdid}} @{$tmpcrd{$crdid}}[-1] \n";
	}
	$pertavg=$perttt/$ncrd;
	my $gr=&gyrationradio(%tmpcrd);
	push (@{$clurank{$cid}},$pertavg,$gr); push(@pertavgall,$pertavg);
    }
    my $rankid=0; my %flag_rank=();undef %flag_rank;
    foreach (sort{$b<=>$a} @pertavgall){
	$rankid++;
	foreach $cid (keys %clurank){
	    if($flag_rank{$cid} ne 1 and abs(@{$clurank{$cid}}[0]-$_) le 0.0001){
		push (@{$clurank{$cid}},$rankid);
		$flag_rank{$cid}=1;
	    }
	}
    }
#    CLURANK: foreach $cid <== $pertavg,$gr,$rankid
#    foreach my $cid (sort {$a<=>$b} keys %clurank) {
#	print "ClusterRANK> $cid -- @{$clurank{$cid}}\n";
#    } exit();
#
    return \%clurank;
    }
sub gyrationradio{
    my %crd=@_;
    my $gr=-1; my @center=(); my $ntt=0;
    foreach my $i (sort {$a<=>$b}keys %crd){
	$ntt++;
	for (my $k=0; $k<3; $k++){ @center[$k]+=@{$crd{$i}}[$k];}
	#print "GYARADIO> $i --  @{$crd{$i}}[0] @{$crd{$i}}[1] @{$crd{$i}}[2] \n";
    }
    for (my $k=0; $k<3; $k++){ @center[$k]=@center[$k]/$ntt}
    
    my $tmpv=0;
    foreach my $i (keys %crd){
	for (my $k=0; $k<3; $k++){ $tmpv+=(@{$crd{$i}}[$k]-@center[$k])**2;}}

    $tmpv=$tmpv/$ntt;
    $gr=sqrt($tmpv); $density=9999;
    $density=(3/4./3.14159)*$ntt/($gr**3) if $gr > 0;
    
    #print "GYARADIO> $gr, $density  -- $ntt\n";

    return $gr;
}





sub parse_arg{
    $numArgs = $#ARGV + 1;
    $index=-1;
    while (++$index<$numArgs){
	$para=$ARGV[$index];
	if($para eq "-c"){
	    $computer=$ARGV[++$index];} 
	elsif($para eq "-p"){
	    $pdblist='-p';
	    $pdbid=$ARGV[++$index];} 
	elsif($para eq "-f"){
	    $pdblist='-f';
	    $pdblistf=$ARGV[++$index];
	    $liststart=$ARGV[++$index];
	    $liststep=$ARGV[++$index];}
	elsif($para eq "-lgsize"){
	    $lgatomnocutoff=$ARGV[++$index];
	    $lgcntatomnocutoff=$ARGV[++$index];}    #ignor ligands with smaller atom numbers
	elsif($para eq "-wcpdb"){
	    $flag_wclusterpdb=1;}
	elsif($para eq "-walldpa"){
	    $flag_walldpa=1;}
	elsif($para eq "-showfit"){
	    $flag_showfit='show';}
	elsif($para eq "-nsec"){
	    $ndpasection=$ARGV[++$index];}
	elsif($para eq "-topp"){
	    $toppercent=$ARGV[++$index];
	    if($toppercent > 1 or $toppercent < 0.1){
		print "Error input for -topp \n";
		exit();}}
	elsif($para eq "-lgcutoff"){
	    $lgbindcutoff=$ARGV[++$index];}
	elsif($para eq "-dpacutoff"){
	    $dpabindingcutoff=$ARGV[++$index];}
	elsif($para eq "-cutoff"){
	    $lgbindcutoff=$ARGV[++$index];
	    $dpabindingcutoff=$ARGV[++$index];}
	elsif($para eq "-clu"){
	    $epislon_dpa=$ARGV[++$index];
	    $MinPts_dpa=$ARGV[++$index]; 
	    $epislon_cluster_dpa=$ARGV[++$index];} 
	elsif($para eq "-het"){
	    $hetatm_name=$ARGV[++$index];}
	elsif($para eq "-adjtopp"){
	    $adjcutperct="true";}
        elsif($para eq "-nolig"){
            $nolig="true";}
	else{
	    print "bad parameter for $para\n";
	    print "$numArgs -- $index : $pdbidlistf,$DATASTART,$DATASTEP,$pdbid,$generateDPA,$DPASTEP,$delta_cutcc,$ev_threshold0,$dpadatadir \n";
	    exit();}
    }
#
#--    
    @revisedRES=('ACE','CGU','MSE','CME','CSS','KCX','TRO','SEP',' FE','OXY'); #%%%%%%
    $hetatm_name='***' if $hetatm_name eq '';  #All the HETATM as LIGAND

    $epislon_dpa=6;
    $MinPts_dpa=3;
    $epislon_cluster_dpa=6;
    $lgatomnocutoff=1;    #ignor ligands with smaller atom numbers
    $lgcntatomnocutoff=1;
    $lgbindcutoff= 6 if ($lgbindcutoff eq '');
    $dpabindingcutoff=6 if ($dpabindingcutoff eq '');

    $toppercent=0.96 if $toppercent eq '';
    $adjcutperct="false" if $adjcutperct eq '';
    $ndpasection=40 if $ndpasection eq '';

    $nolig="false" if $nolig eq '';

    # Updated paths for macOS environment
    if (exists $ENV{DPA_HOME}) {
      $home=$ENV{DPA_HOME};
    } else {
      print "Must set DPA_HOME to repository directory\n";
      exit(1);
    }
    $structuredatadir=$home.'/structuredata';
    $dpadatadir=$home.'/dpadata';
    $exedir=$home.'/bin';
    $workdir=$home.'/scratch';

    $flag_wclusterpdb=-1 if $flag_wclusterpdb eq '';
    $flag_walldpa=-1 if $flag_walldpa eq '';
}
