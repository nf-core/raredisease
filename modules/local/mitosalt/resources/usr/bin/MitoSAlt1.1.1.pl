#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw(first min max sum);
use Data::Dumper;

#INPUT
my $config_file = $ARGV[0];
my $p1 = $ARGV[1];
my $p2 = $ARGV[2];
my $tag = $ARGV[3];

#CREATE DIRECTORIES -- added by ID
mkdir "log" unless -d "log";
mkdir "indel" unless -d "indel";
mkdir "bam" unless -d "bam";
mkdir "tab" unless -d "tab";
mkdir "bw" unless -d "bw";
mkdir "plot" unless -d "plot";

#LOG
open (STDOUT, "| tee -ai log/$tag.log");

my $usage = "usage: perl mitopoint.pl <config_file> <fastq file 1> <fastq file 2> <study name>";
die "Configuration file missing, $usage" unless $config_file;
die "Pair 1 file missing, $usage" unless $p1;
die "Pair 2 file missing, $usage" unless $p2;
die "Study name not given, $usage" unless $tag;

die "check path to fastq 1, script exit.\n" unless -e $p1;
die "check path to fastq 2, script exit.\n" unless -e $p2;


#LOAD CONFIGURATION FILE
open(CONFIG, "<$config_file");
my $User_Preferences;
while (<CONFIG>) {
    chomp;                  # no newline
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences->{$var} = $value;
    my $line = $_;
    print $line."\n";
}

#PROGRAMS
my $hisat2 = $User_Preferences->{hisat2};
my $lastal = $User_Preferences->{lastal};
my $lastsp = $User_Preferences->{lastsp};
my $mfcv = $User_Preferences->{mfcv};
my $b2fq = $User_Preferences->{b2fq};
my $reformat = $User_Preferences->{reformat};
my $samtools = $User_Preferences->{samtools};
my $sambamba = $User_Preferences->{sambamba};
my $gcov = $User_Preferences->{gcov};
my $intersectBed = $User_Preferences->{intersectBed};
my $sortBed = $User_Preferences->{sortBed};
my $clusterBed = $User_Preferences->{clusterBed};
my $randomBed = $User_Preferences->{randomBed};
my $groupBy = $User_Preferences->{groupBy};
my $bg2bw = $User_Preferences->{bg2bw};

#DATABASES
my $hsindex = $User_Preferences->{hsindex};
my $faindex = $User_Preferences->{faindex};
my $lastindex = $User_Preferences->{lastindex};
my $mtfaindex = $User_Preferences->{mtfaindex};
my $gsize = $User_Preferences->{gsize};
my $MT_fasta = $User_Preferences->{MT_fasta};

#COMPUTATION
my $threads = $User_Preferences->{threads};

#MITOCHONDRIA FEATURES
my $refchr = $User_Preferences->{refchr};
my $msize = $User_Preferences->{msize};
my $orihs = $User_Preferences->{orihs};
my $orihe = $User_Preferences->{orihe};
my $orils = $User_Preferences->{orils};
my $orile = $User_Preferences->{orile};
my $exclude = $User_Preferences->{exclude};
my $dloop1 = $exclude;
my $dloop2 = $msize - $exclude;
my $hash;

#SCORING AND FILTERING FEATURES
my $score_threshold = $User_Preferences->{score_threshold};
my $evalue_threshold = $User_Preferences->{evalue_threshold};
my $split_length = $User_Preferences->{split_length};
my $paired_distance = $User_Preferences->{paired_distance};
my $deletion_threshold_min = $User_Preferences->{deletion_threshold_min};
my $deletion_threshold_max = $User_Preferences->{deletion_threshold_max};
my $breakthreshold = $User_Preferences->{breakthreshold};
my $cluster_threshold =  $User_Preferences->{cluster_threshold};
my $breakspan = $User_Preferences->{breakspan};
my $sizelimit = $User_Preferences->{sizelimit};
my $hplimit = $User_Preferences->{hplimit};
my $flank = $User_Preferences->{flank};
my $split_distance_threshold = $User_Preferences->{split_distance_threshold};

#STEPS
my $dna = $User_Preferences->{dna}; #USE SPECIFIC HISAT PARAMETERS FOR DNA AND RNA SEQUENCING
my $enriched = $User_Preferences->{enriched}; #IF THE SEQUENCING IS MITOCHONDRIAL DNA ENRICHED THEN SKIP THE INITIAL NUCLEAR GENOME ALIGNMENT STEP
my $nu_mt = $User_Preferences->{nu_mt}; #MAPPING TO NUCLEAR AND MITOCHONDRIAL GENOME WITH BOWTIE2
my $rmtmp = $User_Preferences->{rmtmp}; #REMOVE TEMPORARY FILES
my $o_mt = $User_Preferences->{o_mt}; #MITOCHONDRIAL READ EXTRACTION AND REMAPPING TO ONLY MITOCHONDRIAL GENOME WITH LASTAL
my $i_del = $User_Preferences->{i_del}; #IDENTIFICATION AND CLUSTERING OF DELETIONS/DUPLICATIONS
my $cn_mt = $User_Preferences->{cn_mt}; #ESTIMATION OF MT DNA COPY NUMBER

#EXIT IF THE CONFIGURATION STEPS DO NOT MATCH
if(($nu_mt eq 'yes' && $enriched eq 'yes')||($nu_mt eq 'no' && $o_mt eq 'no' && $i_del eq 'no')||($cn_mt eq 'yes' && $nu_mt eq 'no')||($cn_mt eq 'yes' && $o_mt eq 'no')||($i_del eq 'yes' && $o_mt eq 'no' && $nu_mt eq 'yes')){
  print scalar(localtime).": Check config file for correct STEP combination\n";
  die;
}

#MAP TO NU+MT GENOME
if($nu_mt eq 'yes' && $enriched eq 'no'){
  print scalar(localtime).": Map to NU+MT genome\n";
  if($dna eq 'yes'){system("$hisat2 -p $threads -x $hsindex --no-temp-splicesite --no-spliced-alignment --max-intronlen 5000 -1 $p1 -2 $p2 -S tmp_$tag.sam 2>> log/$tag.log");}
  if($dna eq 'no'){system("$hisat2 -p $threads -x $hsindex -1 $p1 -2 $p2 -S tmp_$tag.sam 2>> log/$tag.log");}

  #EXTRACT READS
  print scalar(localtime).": Extract reads\n";
  system("$samtools view -@ $threads -bt $faindex tmp_$tag.sam|$samtools sort -@ $threads -o tmp_$tag.bam -");
  system("$samtools index tmp_$tag.bam");
  system("$samtools idxstats tmp_$tag.bam > indel/$tag.count.txt");
  system("$samtools view -u -@ $threads -f 12 tmp_$tag.bam > tmp1_$tag.bam"); #flag to extract unmapped reads with mate also unmapped
  system("$samtools view -u -@ $threads tmp_$tag.bam $refchr> tmp2_$tag.bam");
  system("$samtools merge -@ $threads tmp_final_$tag.bam tmp1_$tag.bam tmp2_$tag.bam");

  #EXTRACT FASTQ
  print scalar(localtime).": Create FASTQ\n";
  system("$samtools fastq tmp_final_$tag.bam > tmp_$tag.fq"); #ADDS PAIR /1, /2 TO THE READ NAMES HENCE BBMAP REFORMAT NOT REQUIRED
  #system("$reformat in=tmp.fq out=tmp1.fq overwrite=true -uniquenames -Xmx100g &>> log/$tag.log");
}

if($nu_mt eq 'yes' && $o_mt eq 'yes' && $enriched eq 'no'){
  #REMAP ON MT GENOME
  print scalar(localtime).": Map to MT genome\n";
  system("$lastal -Q1 -e80 -P$threads $lastindex tmp_$tag.fq|$lastsp > tmp_$tag.maf");
  system("$mfcv sam -d tmp_$tag.maf|$samtools view -@ $threads -bt $mtfaindex -|$samtools sort -@ $threads -o bam/$tag.bam -");
  system("$samtools index bam/$tag.bam");
  system("$mfcv tab tmp_$tag.maf > tab/$tag.tab");

  #COMPRESS TAB
  system("gzip -f tab/$tag.tab");

  #GENERATE BIGWIG
  print scalar(localtime).": Generate Bigwig\n";
  system("$gcov -split -ibam bam/$tag.bam -bg|sort -k1,1 -k2,2n > tmp_$tag.bg");
  system("$bg2bw tmp_$tag.bg $mtfaindex bw/$tag.bw");
  #unlink("tmp_$tag.bg");
}

if($nu_mt eq 'no' && $o_mt eq 'yes' && $enriched eq 'yes'){
  #REMAP ON MT GENOME
  print scalar(localtime).": Map to MT genome\n";
  system("$reformat in=$p1 in2=$p2 out=tmp_$tag.fq overwrite=true addslash=t trimreaddescription=t spaceslash=f -Xmx100g"); #2>> log/$tag.log");
  system("$lastal -Q1 -e80 -P$threads $lastindex tmp_$tag.fq|$lastsp > tmp_$tag.maf");
  system("$mfcv sam -d tmp_$tag.maf|$samtools view -@ $threads -bt $mtfaindex -|$samtools sort -@ $threads -o bam/$tag.bam -");
  system("$samtools index bam/$tag.bam");
  system("$mfcv tab tmp_$tag.maf > tab/$tag.tab");

  #COMPRESS TAB
  system("gzip -f tab/$tag.tab");

  #GENERATE BIGWIG
  print scalar(localtime).": Generate Bigwig\n";
  system("$gcov -split -ibam bam/$tag.bam -bg|sort -k1,1 -k2,2n > tmp_$tag.bg");
  system("$bg2bw tmp_$tag.bg $mtfaindex bw/$tag.bw");
  #unlink("tmp_$tag.bg");
}

if($cn_mt eq 'yes' && $enriched eq 'no' && $o_mt eq 'yes' && $nu_mt eq 'yes'){
  system("$sambamba depth window -w $msize bam/$tag.bam|cut -f 1,5 > indel/$tag.cnmt.txt");
  system("$randomBed -l 1000 -n 3000 -g $gsize > tmp_random_$tag.bed");
  system("$sambamba depth region -L tmp_random_$tag.bed -t $threads -o tmp_random_$tag.cov tmp_$tag.bam");
  system("$groupBy -g 5 -c 8 -o mean -i tmp_random_$tag.cov >> indel/$tag.cnmt.txt");
}

#IDENTIFY DELETIONS/DUPLICATIONS
if($i_del eq 'yes'){
  my $infile = "tab/$tag.tab";
  my $bedfile = "indel/$tag.bed";
  my $breakpointfile = "indel/$tag.breakpoint";
  my $clusterfile = "indel/$tag.cluster";

  my $check_paired = &check_paired($infile);
  print scalar(localtime).": Build split read hash\n";
    $hash = &build_hash($infile,$check_paired);
  &remove_duplicates();

  print scalar(localtime).": Generate non-split read BED\n";
  my $nosplitbed_file = &print_bed($hash,$tag,$refchr);

  print scalar(localtime).": Process hash to get best deletion/duplication candidates\n";
  my $delhash = &process_hash($hash,$bedfile,$breakpointfile,$msize,$tag);

  print scalar(localtime).": Build split read clusters\n";
  my $clusterhash = &get_cluster($delhash,$breakthreshold,$tag);

  print scalar(localtime).": Generate and print results\n";
  my $clustercount = &build_results($hash,$clusterhash,$infile,$refchr,$cluster_threshold,$tag);
  &print_result($clustercount,$clusterfile);

  #PLOT DELETIONS
  print scalar(localtime).":Plot deletions/duplications\n";
  system("R CMD BATCH --no-save --no-restore \'--args $msize $orihs $orihe $orils $orile $sizelimit indel/$tag.cluster indel/$tag.breakpoint $tag $hplimit $MT_fasta $flank\' delplot.R");
}

if($rmtmp eq 'yes'){
  if($nu_mt eq 'yes'){system("rm tmp_$tag* tmp1_$tag* tmp2_$tag* tmp_final_$tag*");}
  if($o_mt eq 'yes' && $enriched eq 'yes'){if(glob("tmp_$tag*")){system("rm tmp_$tag*");}}
  if($i_del eq 'yes'){if(glob("tmp_$tag*")){system("rm tmp_$tag*");}}
  if($cn_mt eq 'yes'){system("rm tmp_random_$tag*");}
}

print scalar(localtime).":Finished\n";
###########################################################SUBROUTINES################################################################
sub check_paired{
  my $check_paired = 0;
  my $path = shift;
  my $line = `zcat $path|tail -n 1`;
  my @elements = split(/\t/, $line);

  $check_paired = 1 if $elements[6]=~m/\/.$/;
  return $check_paired;
}

sub build_hash{
  my $file = shift;
  my $check_paired = shift;

  open(IN, "gunzip -c $file |") || die "can’t open pipe to $file";
  my $hash;
  my $check_hash;
  while(<IN>){
    my $line = $_;
    chomp($line);
    next if $line=~m/^\#/;

    my @elements = split(" ",$line);
    my $score = $elements[0];
    my $chr = $elements[1];

    my $start = $elements[2];
    my $length = $elements[3];
    my $end = $start + $length;

    my $read_start = $elements[7];
    my $id;
    my $mate = 1;

    if($check_paired==0){
      $id = $elements[6];
      if($id=~m/\_.$/){
        $mate = $id;
        $mate =~s/.*\_(.)$/$1/;
        $id =~s/\_.$//;
      }
    }
    elsif($check_paired == 1){
      my @name = split(/\//,$elements[6]);       $id = $name[0];
      $mate = $name[1];
    }

    my $strand = $elements[9];
    my $evalue = $elements[12];
    $evalue=~s/.*\=//;
    next if exists $check_hash->{$id}->{$start.$length};
    $check_hash->{$id}->{$start.$length}++;
    $mate = 2 if $mate > 1;
    #print "$id\t$mate\t$start\t$end\t$strand\t$len\t$score\t$eval\n";

    #FILTER BY SCORE AND MAP EVALUE
    next if $score < $score_threshold;
    next if $evalue > $evalue_threshold;

    push(@{$hash->{$id}->{$mate}->{starts}},$start);
    push(@{$hash->{$id}->{$mate}->{ends}},$end);
    push(@{$hash->{$id}->{$mate}->{rstarts}},$read_start);
    push(@{$hash->{$id}->{$mate}->{strands}},$strand);
    push(@{$hash->{$id}->{$mate}->{scores}},$score);
    push(@{$hash->{$id}->{$mate}->{evalue}},$evalue);
    push(@{$hash->{$id}->{$mate}->{lengths}},$length);
  }
  undef $check_hash;
  close(IN);
  return $hash;
}

sub remove_duplicates{
  my $check_duplicates;

  for my $id (keys %{$hash}){
    my $starts1='NULL';
    my $ends1='NULL';
    my $lengths1='NULL';
    my $starts2='NULL';
    my $ends2='NULL';
    my $lengths2='NULL';
    my $count1 = 0;
    my $count2 = 0;

    if(exists $hash->{$id}->{1}){
      $starts1 = join("",@{$hash->{$id}->{1}->{starts}});
      $ends1 = join("",@{$hash->{$id}->{1}->{ends}});
      $lengths1 = join("",@{$hash->{$id}->{1}->{lengths}});
      $count1 = @{$hash->{$id}->{1}->{starts}};

    }
    if(exists $hash->{$id}->{2}){
      $starts2 = join("",@{$hash->{$id}->{2}->{starts}});
      $ends2 = join("",@{$hash->{$id}->{2}->{ends}});
      $lengths2 = join("",@{$hash->{$id}->{2}->{lengths}});
      $count2 = @{$hash->{$id}->{2}->{starts}};
    }

    my $signature1 = $starts1.":".$ends1.":".$lengths1;
    my $signature2 = $starts2.":".$ends2.":".$lengths2;

    my @signatures = ($signature1,$signature2);
    @signatures = sort @signatures;
    my $signature = join("",@signatures);


    if(exists $check_duplicates->{$signature} && $count1>0 && $count2>0){
      delete $hash->{$id};
    }
    else{
      $check_duplicates->{$signature}++;
    }
  }
}

#PRINT NON-SPLIT BED
sub print_bed{
  my $hash = shift;
  my $tag = shift;
  my $refchr = shift;
  my $filename  = "tmp_".$tag."_nosplit.bed";
  open(ABED,">$filename");
  for my $id (keys %{$hash}){
     #print $id."\n";

    for my $read(keys %{$hash->{$id}}){
      my $name = $id."_".$read;

      #CHECK IF READ IS SPLIT ALIGNED
      my $count = @{$hash->{$id}->{$read}->{starts}};
      next if $count > 1;
      my $start = @{$hash->{$id}->{$read}->{starts}}[0];
      my $end = @{$hash->{$id}->{$read}->{ends}}[0];
      print ABED "$refchr\t$start\t$end\t$name\n";
    }
  }
  system("sort -k2,2n $filename -o $filename");
  close(ABED);
  return $filename;
}

sub process_hash{
  my $hash = shift;
  my $bedfile = shift;
  my $breakpointfile = shift;
  my $size = shift;
  my $tag = shift;
  my $bps = "tmp_".$tag."_bps.bed";
  my $bpe = "tmp_".$tag."_bpe.bed";
  my $delhash=&process_hash1($hash,$bedfile,$breakpointfile,$msize,$bps,$bpe);

  for my $id (keys %{$hash}){
    #print $id."\n";
    for my $read(keys %{$hash->{$id}}){
      #CHECK IF READ IS SPLIT ALIGNED
      my $name = $id."_".$read;
      my $count = @{$hash->{$id}->{$read}->{starts}};
      next unless $count == 2;
      #next unless $name eq 'C1LUFACXX130228:4:1312:19820:41898_1';
      #CHECK IF THE SPLIT READ FRAGMENTS ARE ABOVE LENGTH THRESHOLD
      my $min_len = min @{$hash->{$id}->{$read}->{lengths}};
      next unless $min_len >= $split_length;

      #GET DISTANCE BETWEEN THE FRAGMENTS AND THE BREAKPOINTS
      my @read_starts = @{$hash->{$id}->{$read}->{starts}};
      my @read_ends = @{$hash->{$id}->{$read}->{ends}};
      my @read_strands = @{$hash->{$id}->{$read}->{strands}};
      my @read_lengths = @{$hash->{$id}->{$read}->{lengths}};
      my @read_local_starts = @{$hash->{$id}->{$read}->{rstarts}};

      #LOOK FOR SPLIT READS WHICH MAP INVERSELY SPANNING THE D'LOOP
      my $read_check = 'no';
      if($read_starts[0]<$read_starts[1] && $read_local_starts[0]>$read_local_starts[1]){$read_check = 'yes';}
      elsif($read_starts[0]>$read_starts[1] && $read_local_starts[0]<$read_local_starts[1]){$read_check = 'yes';}
      elsif($read_starts[0]<=$read_starts[1] && $read_local_starts[0]<$read_local_starts[1] && $read_starts[1]<$read_ends[0] ){$read_check = 'yes';}
      elsif($read_starts[0]>=$read_starts[1] && $read_local_starts[0]>$read_local_starts[1] && $read_ends[1]>$read_starts[0] ){$read_check = 'yes';}

      #CHECK IF THE SPLIT READ FRAGMENTS SPAN THE DLOOP
      my $min_start = min @{$hash->{$id}->{$read}->{starts}};
      my $max_end = max @{$hash->{$id}->{$read}->{ends}};
      next if $min_start <= $dloop1 && $max_end >= $dloop2;

      my ($size,$start,$end) = &get_frag_distance(\@read_starts,\@read_ends,$read_check,$msize);
      next if $size < $deletion_threshold_min;
      next if $size > $deletion_threshold_max;
      next if $read_strands[0] ne $read_strands[1];#SPLIT READS IN OPPOSITE STRANDS MAY INDICATE INVERSION
      #print "$name\t$start\t$end\t$size\n";

      $start = 1 if $start == 0;
      $end = 1 if $end == 0;

      #CHECK IF PAIRED READ EXISTS AND ITS ALIGNMENT POSITION SUPPORTS THE SPLIT READ
      my $pair = 1 if $read == 2;
      $pair = 2 if $read == 1;
      my $paired_support = 'no';
      my $distance_paired_support = &paired_support($hash,$id,$read,$pair);
      $paired_support = 'yes' if $distance_paired_support <= $paired_distance;

      #LOOK FOR SPLIT READS WITH LARGE UNMAPPED AREA
      my ($split_distance) = &get_split_distance(\@read_local_starts,\@read_lengths);
      next if $split_distance > $split_distance_threshold;

      #PRINT SPLIT READ IN BED FORMAT
      my @read_scores = @{$hash->{$id}->{$read}->{scores}};
      my @len = &generate_bed($id,$read,\@read_starts,\@read_ends,\@read_lengths,\@read_scores,$read_strands[0]);

      my $readid = $id."_".$read;
      my $clusterid = "cluster_".$count;
      $delhash->{$readid}->{breakstart}=$len[2];
      $delhash->{$readid}->{breakend}=$len[3];
      $delhash->{$readid}->{breaksize}=$len[4];
      $delhash->{$readid}->{readcheck}=$read_check;
      $delhash->{$readid}->{lenstart}=$len[0];
      $delhash->{$readid}->{lenend}=$len[1];

      $start = $len[2];
      $end = $len[3];
      $size = $len[4];

      print BP "$refchr\t$name\t$size\t$start\t$end\t$read_lengths[0]\t$read_lengths[1]\t$paired_support\t$distance_paired_support\t$read_check\n";
      print BPS "$refchr\t$start\t$start\t$name\t0\t+\n" if $read_check eq 'no';
      print BPS "$refchr\t$start\t$start\t$name\t0\t-\n" if $read_check eq 'yes';
      print BPE "$refchr\t$end\t$end\t$name\t0\t+\n" if $read_check eq 'no';
      print BPE "$refchr\t$end\t$end\t$name\t0\t-\n" if $read_check eq 'yes';
    }
  }
  close(BP);
  close(BED);
  close(BPS);
  close(BPE);
  return $delhash;
}

sub process_hash1{
  my $hash = shift;
  my $bedfile = shift;
  my $breakpointfile = shift;
  my $size = shift;
  my $bps = shift;
  my $bpe = shift;

  my $delhash;

  open(BP,">$breakpointfile");
  open(BED,">$bedfile");
  open(BPS,">$bps");
  open(BPE,">$bpe");

  for my $id (keys %{$hash}){
    #print $id."\n";
    for my $read(keys %{$hash->{$id}}){
      #CHECK IF READ IS SPLIT ALIGNED
      my $name = $id."_".$read;
      my $count = @{$hash->{$id}->{$read}->{starts}};
      next unless $count == 3;

      #GET DISTANCE BETWEEN THE FRAGMENTS AND THE BREAKPOINTS
      my @read_starts = @{$hash->{$id}->{$read}->{starts}};
      my @read_ends = @{$hash->{$id}->{$read}->{ends}};
      my @read_strands = @{$hash->{$id}->{$read}->{strands}};
      my @read_lengths = @{$hash->{$id}->{$read}->{lengths}};
      my @read_scores = @{$hash->{$id}->{$read}->{scores}};
      my @read_local_starts = @{$hash->{$id}->{$read}->{rstarts}};


      my %strands=map {$_ => 1} @read_strands;
      my $count_read_strands = keys(%strands);
      next unless $count_read_strands == 1;

      my $as=join(';',@read_starts);my $ae=join(';',@read_ends);my $asr=join(';',@read_strands);my $als=join(';',@read_local_starts);
      my @pos = (1000,1000,1000,1000);
      #CHECK FOR FALSE SPLITS AT ORIGIN WITH A REAL DELETION
      $pos[0] = 0 if grep {$_ == 0} @read_starts;
      $pos[1] = $msize if grep {$_ == $msize} @read_starts;
      $pos[2] = 0 if grep {$_ == 0} @read_ends;
      $pos[3] = $msize if grep {$_ == $msize} @read_ends;

      #print "$name\t$as\t$ae\t$asr\t$als\t$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\n";
      if(($pos[0] != 1000||$pos[1] != 1000)&&($pos[2] != 1000||$pos[3] != 1000)){
        my $sort_hash;
        for(my $i=0;$i<3;$i++){
          $sort_hash->{$read_local_starts[$i]}->{start}=$read_starts[$i];
          $sort_hash->{$read_local_starts[$i]}->{end}=$read_ends[$i];
          $sort_hash->{$read_local_starts[$i]}->{strand}=$read_strands[$i];
          $sort_hash->{$read_local_starts[$i]}->{length}=$read_lengths[$i];
          $sort_hash->{$read_local_starts[$i]}->{score}=$read_scores[$i];
        }
        my $order_count=1;
        foreach my $local_start (sort { $a <=> $b } keys %{$sort_hash}){
          my $start_check=$sort_hash->{$local_start}->{start};
          my $end_check=$sort_hash->{$local_start}->{end};
          my $pos_check='no';
          $pos_check='yes' unless($start_check==0||$start_check==$msize||$end_check==0||$end_check==$msize);
          if($pos_check eq 'yes'){last;}
          $order_count++;
        }

        my $delete_count=3 if $order_count==1;
        $delete_count=1 if $order_count==3;
        $order_count=1;
        foreach my $local_start (sort { $a <=> $b } keys %{$sort_hash}){
          delete $sort_hash->{$local_start} if($order_count==$delete_count);
          $order_count++;
        }
        my @read_startsN;
        my @read_endsN;
        my @read_strandsN;
        my @read_local_startsN;
        my @read_scoresN;
        my @read_lengthsN;
        foreach my $local_start (sort { $a <=> $b } keys %{$sort_hash}){
          my $start =   $sort_hash->{$local_start}->{start};
          my $end =   $sort_hash->{$local_start}->{end};
          my $score =   $sort_hash->{$local_start}->{score};
          my $length =   $sort_hash->{$local_start}->{length};
          my $strand =   $sort_hash->{$local_start}->{strand};
          push @read_startsN,$start;
          push @read_endsN,$end;
          push @read_local_startsN,$local_start;
          push @read_strandsN,$strand;
          push @read_scoresN,$score;
          push @read_lengthsN,$length;
        }

        #CHECK IF THE SPLIT READ FRAGMENTS ARE ABOVE LENGTH THRESHOLD
        my $min_len = min @read_lengthsN;
        next unless $min_len >= $split_length;


        #LOOK FOR SPLIT READS WHICH MAP INVERSELY SPANNING THE D'LOOP
        my $read_check = 'no';
        if($read_startsN[0]<$read_startsN[1] && $read_local_startsN[0]>$read_local_startsN[1]){$read_check = 'yes';}
        elsif($read_startsN[0]>$read_startsN[1] && $read_local_startsN[0]<$read_local_startsN[1]){$read_check = 'yes';}
        my ($size,$start,$end) = &get_frag_distance(\@read_startsN,\@read_endsN,$read_check,$msize);
        next if $size < $deletion_threshold_min;
        next if $size > $deletion_threshold_max;

        #CHECK IF PAIRED READ EXISTS AND ITS ALIGNMENT POSITION SUPPORTS THE SPLIT READ
        my $pair = 1 if $read == 2;
        $pair = 2 if $read == 1;
        my $paired_support = 'no';
        my $distance_paired_support = &paired_support($hash,$id,$read,$pair);
        $paired_support = 'yes' if $distance_paired_support <= $paired_distance;

        $start = 1 if $start == 0;
        $end = 1 if $end == 0;

        #LOOK FOR SPLIT READS WITH LARGE UNMAPPED AREA
        my ($split_distance) = &get_split_distance(\@read_local_startsN,\@read_lengthsN);
        next if $split_distance > 5;

        #PRINT SPLIT READ IN BED FORMAT
        my @len = &generate_bed1($id,$read,\@read_startsN,\@read_endsN,\@read_lengthsN,\@read_scoresN,$read_strandsN[0]);
        my $readid = $id."_".$read;
        my $clusterid = "cluster_".$count;
        $delhash->{$readid}->{breakstart}=$start;
        $delhash->{$readid}->{breakend}=$end;
        $delhash->{$readid}->{breaksize}=$size;
        $delhash->{$readid}->{readcheck}=$read_check;
        $delhash->{$readid}->{lenstart}=$len[0];
        $delhash->{$readid}->{lenend}=$len[1];

        #$start = $len[2];
        #$end = $len[3];
        #$size = $len[4];

        print BP "$refchr\t$name\t$size\t$start\t$end\t$read_lengths[0]\t$read_lengths[1]\t$paired_support\t$distance_paired_support\t$read_check\n";
        print BPS "$refchr\t$start\t$start\t$name\t0\t+\n" if $read_check eq 'no';
        print BPS "$refchr\t$start\t$start\t$name\t0\t-\n" if $read_check eq 'yes';
        print BPE "$refchr\t$end\t$end\t$name\t0\t+\n" if $read_check eq 'no';
        print BPE "$refchr\t$end\t$end\t$name\t0\t-\n" if $read_check eq 'yes';

      }
    }
  }
  return $delhash;
}

#CHECK IF THE PAIR OF A SPLIT READ LIES WITHIN A GIVEN THRESHOLD DISTANCE TO A FRAGMENT OF THE SPLIT READ
sub paired_support{
   my $hash = shift;
   my $id = shift;
   my $read = shift;
   my $pair = shift;
   my $distance = 1000;

   if(exists $hash->{$id}->{$pair}){
     my $distance1 = 0;
     my $distance2=0;

     my @read_starts = @{$hash->{$id}->{$read}->{starts}};
     my @read_ends = @{$hash->{$id}->{$read}->{ends}};

     my $pair_start = @{$hash->{$id}->{$pair}->{starts}}[0];
     my $pair_end = @{$hash->{$id}->{$pair}->{ends}}[0];

     $distance1 = $read_starts[0] - $pair_end if $read_starts[0] > $pair_end;
     $distance1 = $pair_start - $read_ends[0] if $pair_start > $read_ends[0];

     $distance2 = $read_starts[1] - $pair_end if $read_starts[1] > $pair_end;
     $distance2 = $pair_start - $read_ends[1] if $pair_start > $read_ends[1];
     $distance = min ($distance1,$distance2);
   }
   return $distance;
}

#GET THE SIZE OF PUTATIVE DELETION AND THE BREAKPOINTS
sub get_frag_distance{
  my $starts = shift;
  my $ends = shift;
  my $read_check = shift;
  my $msize = shift;
  my @res;

  my $frag1_start = @{$starts}[0];
  my $frag2_start = @{$starts}[1];

  my $frag1_end = @{$ends}[0];
  my $frag2_end = @{$ends}[1];

  if($read_check eq 'no'){
    if($frag1_start > $frag2_end){
      my $size = $frag1_start - $frag2_end;
      @res = ($size,$frag2_end,$frag1_start);
    }
    elsif($frag2_start > $frag1_end){
      my $size = $frag2_start - $frag1_end;
      @res = ($size,$frag1_end,$frag2_start);
    }
    elsif($frag1_start < $frag2_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag2_start < $frag1_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag1_start <= $frag2_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_start + $frag1_end ;
      @res = ($size,$frag1_end,$frag2_start);
    }
    elsif($frag2_start <= $frag1_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag1_end + $frag2_start ;
      @res = ($size,$frag2_start,$frag1_end);
    }
  }

  if($read_check eq 'yes'){
    if($frag1_start > $frag2_end){
      my $size = ($msize -$frag1_end) + $frag2_start;
      @res = ($size,$frag2_start,$frag1_end);
    }
    elsif($frag2_start > $frag1_end){
      my $size = ($msize-$frag2_end) + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag1_start <= $frag2_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_start + $frag1_end ;
      @res = ($size,$frag1_end,$frag2_start);
    }
    elsif($frag2_start <= $frag1_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag1_end + $frag2_start ;
      @res = ($size,$frag2_start,$frag1_end);
    }
    elsif($frag1_start < $frag2_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag2_start < $frag1_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
  }
  return @res;
}

#GET THE DISTANCE BETWEEN SPLIT READS
sub get_split_distance{
  my $starts = shift;
  my $lengths = shift;
  my @res;

  my $local1_start = @{$starts}[0];
  my $local2_start = @{$starts}[1];

  my $local1_length = @{$lengths}[0];
  my $local2_length = @{$lengths}[1];

  if($local1_start>$local2_start){
    my $distance = $local1_start - $local2_start - $local2_length;
    @res = ($distance);
  }
  if($local2_start>$local1_start){
    my $distance = $local2_start - $local1_start - $local1_length;
    @res = ($distance);
  }
  return @res;
}

#GENERATE A BED FILE FOR THE SPLIT READS FOR IGV VISUALIZATION
sub generate_bed{
  my $id = shift;
  my $read = shift;
  my $starts = shift;
  my $ends = shift;
  my $lengths = shift;
  my $scores = shift;
  my $strand = shift;

  my $size = shift;
  my $bstart = shift;
  my $bend = shift;
  my @len;

  my $start;
  my $end;

  my $frag1_start = @{$starts}[0];
  my $frag2_start = @{$starts}[1];

  my $frag1_end = @{$ends}[0];
  my $frag2_end = @{$ends}[1];

  my $frag1_length = @{$lengths}[0];
  my $frag2_length = @{$lengths}[1];

  my $frag1_score = @{$scores}[0];
  my $frag2_score = @{$scores}[1];
  my $score = int(($frag1_score+$frag2_score)/2);

  if($frag1_start > $frag2_end && $strand eq '-'){
    $start = $frag2_start;
    $end = $frag1_end;
    $bstart = $frag2_end;
    $bend = $frag1_start;
    $size = $frag1_start - $frag2_end;
    my $block_start=$frag1_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag2_length\,$frag1_length\t0\,$block_start\n";
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }
  elsif($frag1_start > $frag2_end && $strand eq '+'){
    $start = $frag2_start;
    $end = $frag1_end;
    $bstart = $frag2_start;
    $bend = $frag1_end;
    $size = 1 + $msize- $frag1_end + $frag2_start;
    my $block_start=$frag1_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag2_length\,$frag1_length\t0\,$block_start\n";
    @len = ($frag1_length,$frag2_length,$bstart,$bend,$size);
  }
  elsif($frag2_start > $frag1_end && $strand eq '+'){
    $start = $frag1_start;
    $end = $frag2_end;
    $bstart = $frag1_end;
    $bend = $frag2_start;
    $size = $frag2_start - $frag1_end;
    my $block_start=$frag2_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";
    @len = ($frag1_length,$frag2_length,$bstart,$bend,$size);
  }
  elsif($frag2_start > $frag1_end && $strand eq '-'){
    $start = $frag1_start;
    $end = $frag2_end;
    $bstart = $frag1_start;
    $bend = $frag2_end;
    $size = 1 + $msize - $frag2_end + $frag1_start;
    my $block_start=$frag2_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }
  elsif($frag2_start < $frag1_end && $strand eq '+'){
    $start = $frag1_start;
    $end = $frag2_end;
    $bstart = $frag2_start;
    $bend = $frag1_end;
    $size = $msize - $frag1_end + $frag2_start;
    my $block_start=$frag1_end-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }
  elsif($frag1_start < $frag2_end && $strand eq '-'){
    $start = $frag2_start;
    $end = $frag1_end;
    $bstart = $frag1_start;
    $bend = $frag2_end;
    $size = $msize - $frag2_end + $frag1_start;
    my $block_start=$frag2_end-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }

  return @len;
}

#GENERATE A BED FILE FOR THE SPLIT READS FOR IGV VISUALIZATION (reads which split twice, once due to deletion/duplication and once due to genome circularity)
sub generate_bed1{
  my $id = shift;
  my $read = shift;
  my $starts = shift;
  my $ends = shift;
  my $lengths = shift;
  my $scores = shift;
  my $strand = shift;
  my @len;

  my $start;
  my $end;

  my $frag1_start = @{$starts}[0];
  my $frag2_start = @{$starts}[1];

  my $frag1_end = @{$ends}[0];
  my $frag2_end = @{$ends}[1];

  my $frag1_length = @{$lengths}[0];
  my $frag2_length = @{$lengths}[1];

  my $frag1_score = @{$scores}[0];
  my $frag2_score = @{$scores}[1];
  my $score = int(($frag1_score+$frag2_score)/2);

  if($frag1_start > $frag2_end){
    $start = $frag2_start;
    $end = $frag1_end;
    my $block_start=$frag1_start-$start;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag2_length\,$frag1_length\t0\,$block_start\n";
    @len = ($frag2_length,$frag1_length);
  }
  elsif($frag2_start > $frag1_end){
    $start = $frag1_start;
    $end = $frag2_end;
    my $block_start=$frag2_start-$start;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";
    @len = ($frag1_length,$frag2_length);
  }
  return @len;
}

sub get_cluster{
  my $delhash = shift;
  my $breakthreshold = shift;
  my $tag = shift;
  my $bps = "tmp_".$tag."_bps.bed";
  my $bpe = "tmp_".$tag."_bpe.bed";
  my $bpsc = "tmp_".$tag."_bps.cls";
  my $bpec = "tmp_".$tag."_bpe.cls";
  my $clusterhash;

  system("$sortBed -i $bps|$clusterBed -s -d $breakthreshold -i stdin > $bpsc");
  system("$sortBed -i $bpe|$clusterBed -s -d $breakthreshold -i stdin > $bpec");

  open(BPS, "<$bpsc");
  my @bps=<BPS>;
  my $bps_hash;
  foreach my $line(@bps){
    chomp $line;
    my @fields=split(/\t/,$line);
    my $readid=$fields[3];
    my $clusterid=$fields[6];
    $bps_hash->{$readid} = $clusterid;
  }
  close(BPS);

  open(BPE, "<$bpec");
  my @bpe=<BPE>;
  my $bpe_hash;
  foreach my $line(@bpe){
    chomp $line;
    my @fields=split(/\t/,$line);
    my $readid=$fields[3];
    my $clusterid=$fields[6];
    $bpe_hash->{$readid} = $clusterid;
  }
  close(BPE);

  for my $readid(keys %{$delhash}){
    next unless exists $bps_hash->{$readid};
    next unless exists $bpe_hash->{$readid};
    my $clusterids = $bps_hash->{$readid};
    my $clusteride = $bpe_hash->{$readid};
    my $clusterid="cluster".$clusterids.$clusteride;
    $clusterhash->{$clusterid}->{$readid} = $delhash->{$readid};
  }
  return $clusterhash;
}

sub build_results{
  my $hash = shift;
  my $clusterhash = shift;
  my $infile = shift;
  my $refchr = shift;
  my $cluster_read_count_threshold = shift;
  my $tag = shift;
  my $clustercount;

  my $nosplitbed_file = "tmp_".$tag."_nosplit.bed";;
  for my $clusterid (keys %{$clusterhash}){
    my $cluster_read_count = scalar keys %{$clusterhash->{$clusterid}};
    next if $cluster_read_count < $cluster_read_count_threshold;

    print "Check $clusterid\n";
    my $clustercheck;
    $clustercount->{$clusterid}->{wt} = 0;
    $clustercount->{$clusterid}->{mt} = scalar keys %{$clusterhash->{$clusterid}};

    my $splitbed_fileS = "tmp_".$tag.".split.start.bed";
    my $splitbed_fileE = "tmp_".$tag.".split.end.bed";

    open(SBEDS,">$splitbed_fileS");
    open(SBEDE,">$splitbed_fileE");

    for my $readid (keys %{$clusterhash->{$clusterid}}){
      my $breakstart = $clusterhash->{$clusterid}->{$readid}->{breakstart};
      my $breakend = $clusterhash->{$clusterid}->{$readid}->{breakend};
      my $lenstart = $clusterhash->{$clusterid}->{$readid}->{lenstart};
      my $lenend = $clusterhash->{$clusterid}->{$readid}->{lenend};

      push(@{$clustercount->{$clusterid}->{starts}},$breakstart);
      push(@{$clustercount->{$clusterid}->{ends}},$breakend);
      push(@{$clustercount->{$clusterid}->{names}},$readid);
      push(@{$clustercount->{$clusterid}->{lenstarts}},$lenstart);
      push(@{$clustercount->{$clusterid}->{lenends}},$lenend);

      print SBEDS "$refchr\t$breakstart\t$breakstart\n";
      print SBEDE "$refchr\t$breakend\t$breakend\n";
    }

    system("sort -u -k2,2n -k3 $splitbed_fileS -o $splitbed_fileS");
    system("sort -u -k2,2n -k3 $splitbed_fileE -o $splitbed_fileE");

    my $intersectbed_fileS = "tmp_".$tag.".intersect.start.bed";
    my $intersectbed_fileE = "tmp_".$tag.".intersect.end.bed";


    system("$intersectBed -wo -sorted -a $nosplitbed_file -b $splitbed_fileS > $intersectbed_fileS");
    open(IS,"<$intersectbed_fileS");
    while(<IS>){
      my $line = $_;
      chomp($line);
      my @crd = split(/\t/,$line);
      my $start = $crd[1];
      my $end = $crd[2];
      my $name = $crd[3];
      my $breakstart = $crd[5];
      next if exists $clustercheck->{$name};
      my $diffstart = $breakstart - $start;
      my $diffend = $end - $breakstart;
      next if $diffstart <= $breakspan || $diffend <= $breakspan;
      $clustercount->{$clusterid}->{wt1}++;
      $clustercheck->{$name}++;
    }
    close(SBEDS);
    close(IS);
    unlink $splitbed_fileS,$intersectbed_fileS;

    system("$intersectBed -wo -sorted -a $nosplitbed_file -b $splitbed_fileE > $intersectbed_fileE");
    open(IE,"<$intersectbed_fileE");
    while(<IE>){
      my $line = $_;
      chomp($line);
      my @crd = split(/\t/,$line);
      my $start = $crd[1];
      my $end = $crd[2];
      my $name = $crd[3];
      my $breakend = $crd[5];
      next if exists $clustercheck->{$name};
      my $diffstart = $breakend - $start;
      my $diffend = $end - $breakend;
      next if $diffstart <= $breakspan || $diffend <= $breakspan;
      $clustercount->{$clusterid}->{wt2}++;
      $clustercheck->{$name}++;
    }
    close(SBEDE);
    close(IE);
    unlink $splitbed_fileE,$intersectbed_fileE;
  }
  return $clustercount;
}

sub print_result{
  my $clustercount = shift;
  my $clusterfile = shift;
  open(CF,">$clusterfile");

  for my $clusterid(keys %{$clustercount}){
    my $names = $clustercount->{$clusterid}->{names};
    my $starts = $clustercount->{$clusterid}->{starts};
    my $ends = $clustercount->{$clusterid}->{ends};
    my $lenstarts = $clustercount->{$clusterid}->{lenstarts};
    my $lenends = $clustercount->{$clusterid}->{lenends};
    my $mt = $clustercount->{$clusterid}->{mt};
    my $wt1 = 0;
    $wt1 = $clustercount->{$clusterid}->{wt1} if exists $clustercount->{$clusterid}->{wt1};
    my $wt2 = 0;
    $wt2 = $clustercount->{$clusterid}->{wt2} if exists $clustercount->{$clusterid}->{wt2};
    my $sum_wt = $wt1 + $wt2;
    my $mean_wt = sprintf "%.0f",$sum_wt/2;

    my $names_print = join("\,",@$names);
    my $starts_print = join("\,",@$starts);
    my $ends_print = join("\,",@$ends);

    #if fragments overlap, then lenstart and and lenend are committed to NA
    my $lenstarts_print = "NA";
    my $lenends_print = "NA";
    if(defined(@$lenstarts[0])){$lenstarts_print = join("\,",@$lenstarts);}
    if(defined(@$lenends[0])){$lenends_print = join("\,",@$lenends);}

    my $perc_hp = $mt*100/($mt+$mean_wt);

    print CF "$clusterid\t$names_print\t$starts_print\t$ends_print\t$lenstarts_print\t$lenends_print\t$mt\t$mean_wt\t$perc_hp\n";
  }
  close(CF);
}
