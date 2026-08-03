#!/usr/bin/perl
#
#-- Based on Backbone-enhanced ENM
#-- Generating DPA data for PDB protein structure
#     To predict Protein/Ligand interaction

#
#-- TOP output file
$proccf='_calc_dpa.rcd';
&parse_arg();

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

foreach $pdbid ( @idlist){

    my ($calc_info,$infocalc,$nca)=&generate_dpa($pdbid,$structuredatadir,$exedir,$workdir,$dpadir,
			$dpastep,$cutcc,$cutlg,$cutll,$wcc,$wctc,$wlg,$wll,$ev_threshold);
    
    if($calc_info !=0) {
	my $timetitle=&timetitle();open(pp,">> ./$proccf"); print pp "$pdbid $infocalc\n";close(pp);}
    print "$pdbid -- $infocalc\n" if($calc_info != 0);
    printf "%4s\t%5d\t%9.4f%12s\n",$pdbid,$nca,$infocalc," (CPU seconds)" if $calc_info==0;
}
exit();


sub generate_dpa{
    my ($pdbid,$structuredatadir,$exedir,$workdir,$dpadir,
	$dpastep,$cutcc,$cutlg,$cutll,$wcc,$wctc,$wlg,$wll,$ev_threshold)=@_;

    my $pdbfile=$pdbid.'.pdb';my $cacrdf=$pdbid.'_ca.crd';my $selechain='ALL';
    
    return -9999,"PDB FILE $pdbfile does not exist in $structuredatadir" if(! -e "$structuredatadir/$pdbfile");
    my ($nca,$tmp)=&read_pdb($structuredatadir,$pdbfile,$selechain,'getcacrd',$cacrdf,$workdir);
    return -1000,' CA_CRD data error' if($nca<=5);
    my $msmssurff=$pdbid.'.surf'; my $dpaf=$pdbid.'.dpa';
    return -9998,"SURF FILE $msmssurff does not exist in $structuredatadir" if(! -e "$structuredatadir/$msmssurff");

    my $cutcc0=$cutcc;my $cutlg0=$cutlg;my $cutll0=$cutll;
    my $wcc0=$wcc;my $wctc0=$wctc;my $wlg0=$wlg; my $wll0=$wll;

    my $dpainputf=$pdbid.'_dpa.inp';my $dpaoutputf=$pdbid.'_dpa.out';
dpacycle:
    open(dpain,">$workdir/$dpainputf");
    print dpain "\"$workdir/$cacrdf\" \"$structuredatadir/$msmssurff\" \"$dpadir/$dpaf\"\n";
    print dpain "$wcc0 \t$wctc0 \t$wlg0 \t$wll0\n";
    print dpain "$cutcc0 \t$cutlg0 \t$cutll0\n";
    print dpain "$ev_threshold\n";
    close(dpain);
    system "$exedir/ggsspnma < $workdir/$dpainputf > $workdir/$dpaoutputf";

#    system "$exedir/gsspnma.exe < $workdir/$dpainputf > $workdir/$dpaoutputf";

    my ($info,$info2)=&read_calc_info($workdir,$dpaoutputf);  #info2: CPU time
    &dpainfo_procc($pdbid,$nca,$info,$info2,".");
    if($info==-1){$cutcc0+=1;$cutlg0=$cutcc0+5; goto dpacycle;}

    system "rm -f $workdir/$dpainputf $workdir/$dpaoutputf $workdir/$cacrdf" if $flag_check eq 'false';

    return 0,$info2,$nca;
}
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
sub parse_arg{
    $ligandatom_numcut=4;    #ignor ligands with smaller atom numbers
    $ligandcontactatom_numcut=4;
    $domainsize_cutoff=1000;
    $lgbinding_cutoff=6;
    $dpabinding_cutoff=6;

    $dens_msms='';
    $epislon=5;
    $MinPts=3;
    $epislon_cluster=5;
    $computer=1;
    $hetatm_name='***';
    my $densDPASTEP;
    my $gdpaDPASTEP;
    $DPAINI=1;
    $numArgs = $#ARGV + 1;
    $index=-1;
    while (++$index<$numArgs){
	$para=$ARGV[$index];
	if($para eq "-p"){
	    $pdblist='-p';
	    $pdbid=$ARGV[++$index];} 
	elsif($para eq "-f"){
	    $pdblist='-f';
	    $pdblistf=$ARGV[++$index];
	    $liststart=$ARGV[++$index];
	    $liststep=$ARGV[++$index];}
	elsif($para eq "-gdpa"){
	    $dpastart=$ARGV[++$index]; 
	    $dpastep=$ARGV[++$index]; 
	    $delta_cutcc=$ARGV[++$index];
	    $ev_threshold=$ARGV[++$index];
	    $deltacutlg=$ARGV[++$index];}
	elsif($para eq "-wcc"){
	    $wcc=$ARGV[++$index];}
	elsif($para eq "-check"){
	    $flag_check='true';}
	elsif($para eq "-wctc"){
	    $wctc=$ARGV[++$index];}
	elsif($para eq "-wlg"){
	    $wlg=$ARGV[++$index];}
	elsif($para eq "-wll"){
	    $wll=$ARGV[++$index];}
	elsif($para eq "-cutcc"){
	    $cutcc=$ARGV[++$index];}
	elsif($para eq "-dcutcc"){
	    $dcutcc=$ARGV[++$index];}
	elsif($para eq "-cutlg"){
	    $cutlg=$ARGV[++$index];}
	elsif($para eq "-cutll"){
	    $cutll=$ARGV[++$index];}
	else{
	    open(pp,">> $procf");
	    my $timetitle = &timetitle(); 
	    print pp "\nRecording at $timetitle\n";
	    print pp "bad parameter for $para\n";
	    close(pp);
	    exit();}
    }

    $liststart=1 if $liststart eq '';
    $liststep=1  if $liststep  eq '';

    $wcc=1. if $wcc eq '';
    $wctc=36. if $wctc eq '';
    $wlg=12.  if $wlg  eq '';
    $wll=10.  if $wll  eq '';
    $cutcc=10 if $cutcc eq '';
    $dcutcc=1. if $dcutcc eq '';
    $cutlg=15.  if $cutlg eq '';
    $cutll=13. if $cutll eq '';
    $ev_threshold=0.001 if $ev_threshold eq '';
    $cutlg=$cutcc+5. if $cutlg < $cutcc+5.;
    $flag_check='false' if $flag_check eq '';
    if (exists $ENV{DPA_HOME}) {
      $home=$ENV{DPA_HOME}
    } else {
      print "Must set DPA_HOME to repository directory\n";
      exit(1);
    }
    $scrh=$home.'/scratch';
    $workdir=$scrh;
    $outputdir=$scrh;
    $exedir=$home.'/bin';
    $structuredatadir=$home.'/structuredata';
    $dpadir=$home.'/dpadata';
}

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
sub dpainfo_procc{
    my ($pdbid,$nca,$info,$info2,$outputdir)=@_;

    my $timetitle = &timetitle();
    open(pp,">> $outputdir/$proccf");
    if($info==0) {
	printf pp "%4s\t%5d\t%9.4f%31s %14s\n",$pdbid,$nca,$info2," (CPU seconds) -- Recording at ",$timetitle;}
    if($info==-9999){
	print pp "$pdbid OPEN cacrd error   -- Recording at $timetitle\n";}
    elsif($info==-9998){
	print pp "$pdbid READ cacrd error   -- Recording at $timetitle\n";}	
    elsif($info==-9997){
	print pp "$pdbid OPEN surfcrd error   -- Recording at $timetitle\n";}
    elsif($info==-9996){
	print pp "$pdbid READ surfcrd error   -- Recording at $timetitle\n";}
    elsif($info==-5000){
	print pp "$pdbid protein CRD pair-distance < 0.3   -- Recording at $timetitle\n";}
    elsif($info==-5100){
	print pp "$pdbid protein/ligand (surf point) CRD pair-distance < 1   -- Recording at $timetitle\n";}
    elsif($info==-2000){
	print pp "$pdbid Dimension of RP (neighboring contacts) is too smal   -- Recording at $timetitle\n";}
    elsif($info==-1000){
	print pp "$pdbid Dimension of HESSIAN is too smal   -- Recording at $timetitle\n";}
    elsif($info==-9){
	print pp "$pdbid Kessian is singular (NO \"G^t K^-1 G\")    -- Recording at $timetitle\n";}
    elsif($info==-1){
	print pp "$pdbid CUTCC for APOprotein is too small    -- Recording at $timetitle\n";}
    elsif($info='null'){
	print pp "$pdbid\tDPA calc. segment false  -- Recording at $timetitle\n";}
    close(pp);
}
sub read_calc_info{
    my $dir=@_[0]; my $file=@_[1];
    my $info='null'; my $info2='null';
    open(tmpf,"<  $dir/$file");
    while(<tmpf>){
	my @tmp=split /\s+/,$_;
	$info=@tmp[1];  $info=@tmp[0] if @tmp[0] ne '';
	$info2=@tmp[2]; $info2=@tmp[1] if @tmp[0] ne '';
    }
    close(tmpf);
    return $info,$info2;
}
