#!/bin/bash

# This script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# function "build" to compile
function build() {
  cargo build --release
  if [ $? -gt 0 ]; then
        echo "ERROR compiling sepia, did you install Rust(https://www.rust-lang.org/)?";
        exit 1
  fi
}

function print_separator() {
  GREEN_LINE="\033[0;32m -------------------------------------------------------------- \033[0m"
  echo -e "${GREEN_LINE}"
}

# DEFAULTS
THREADS=2

# Get options: -b build, -c clean
while getopts "bcht:" opt; do
  case $opt in
    # BUILD
    b)
      echo " [*] Building"
      build
      ;;
    # CLEAN
    c)
      echo " [*] Cleaning..."
      for FILE in target/release demo_index {test_classify,Ecoli,Efaecium}*.{gz,txt};
      do
        
        if [[ -e "$DIR"/"$FILE" ]]; then
          rm -rf "$DIR"/"$FILE"
        fi
      done
      exit
      ;;
    # -t INT (Threads)
    t)
      echo " - Setting THREADS=$OPTARG"
      THREADS=$OPTARG
      ;;
    h)
      echo "Help:"
      echo "  -b:     Force rebuilding the binary"
      echo "  -c:     Clean test.sh products (binaries and output files)"
      echo "  -t INT: Number of threads for multithreading [default: 2]"
      echo "  -h:     help"
      ;;
    \?)
      echo "Invalid option: -$OPTARG"
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      ;;
  esac
done
 


#compile sepia in release mode if the binary is not already present

if [[ ! -e "$DIR"/target/release/sepia ]]; then
  print_separator
  echo "Sepia binary not found...";
  build
fi

# Build the database with one thread
print_separator;
echo " [1] Building the index";
./target/release/sepia build -k 31 -m 21 -i ./demo_index -r ./ref_demo.txt
if [[ $? -gt 0 ]]; then
  echo "ERROR building sepia index ./demo_index using ./ref_demo.txt";
  exit 1
fi

# Build the database with two threads
print_separator
echo " [2] Building the index (multithreads)"
./target/release/sepia build -k 31 -m 21 -i ./demo_index -r ./ref_demo.txt -p $THREADS -c 2
if [[ $? -gt 0 ]]; then
  echo "ERROR building sepia index ./demo_index using ./ref_demo.txt";
  exit 1
fi

#single end read classification
print_separator
echo " [3] Classifying single end reads"
./target/release/sepia classify -i ./demo_index -q ./1000K_test_fq/Ecoli1_1000K_1.fastq.gz -n test_classify -t $THREADS
if [ $? -gt 0 ]; then
  echo "ERROR classifying reads ./test_data/SRR548019.fastq.gz with ./test_data/phage.bxi"
  exit 1
fi

#paired end read classification
print_separator
echo " [4] Classifying paired-end reads"
./target/release/sepia classify -i ./demo_index -q ./1000K_test_fq/Ecoli1_1000K_*.fastq.gz -n test_classify -t $THREADS
if [ $? -gt 0 ]; then
  echo "ERROR classifying reads ./test_data/SRR548019.fastq.gz with ./test_data/phage.bxi"
  exit 1
fi

#Test the output paired end read classification
print_separator
echo " [5] Testing the output";
declare -a expected=('root;d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli_D	56	0.807	24029')
lastIndex=$((${#expected[@]} - 1))
#echo "$lastIndex .. ${expected[@]}"
for i in $(seq 0 $lastIndex); do
  j=$((i + 1));
  #observed="$(grep coli_D test_classify_summary.txt | cut -f $j -d $'\t')"
  observed="$(grep coli_D test_classify_summary.txt)"
  echo "${expected[$i]} <=> $observed";
  if [ "${expected[$i]}" != "$observed" ]; then
    echo "ERROR: test output was incorrect in field $i.  I expected ${expected[$i]} but got $observed."
    exit 1;
  fi
done
echo " [OK] All fields matched!"

#paired end read batch classification (TODO)
print_separator
echo " [6] Classifying paired-end reads"
./target/release/sepia batch_classify -i ./demo_index -q batch_example.txt -T batch_example -t $THREADS
if [ $? -gt 0 ]; then
  echo "ERROR batch classifying reads batch_example.txt with demo_index"
  exit 1
fi

#Test the output paired end read classification
print_separator
echo " [7] Testing the output batch classification";
declare -a expected1=('root;d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli_D	56	0.807	24029')
declare -a expected2=('root;d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli_D	62	0.775	24338')
declare -a expected3=('root;d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia albertii	1	0.038	501')
#lastIndex=$((${#expected[@]} - 1))
#echo "$lastIndex .. ${expected[@]}"
for i in $(seq 0 $lastIndex); do
  j=$((i + 1));
  #observed="$(grep coli_D Ecoli1_batch_example_summary.txt | cut -f $j -d $'\t')"
  observed="$(grep coli_D Ecoli1_batch_example_summary.txt)"
  echo "${expected1[$i]} <=> $observed";
  if [ "${expected1[$i]}" != "$observed" ]; then
    echo "ERROR: test output was incorrect in field $i.  I expected ${expected1[$i]} but got $observed."
    exit 1;
  fi
done
for i in $(seq 0 $lastIndex); do
  j=$((i + 1));
  #observed="$(grep coli_D Ecoli2_batch_example_summary.txt | cut -f $j -d $'\t')"
  observed="$(grep coli_D Ecoli2_batch_example_summary.txt)"
  echo "${expected2[$i]} <=> $observed";
  if [ "${expected2[$i]}" != "$observed" ]; then
    echo "ERROR: test output was incorrect in field $i.  I expected ${expected2[$i]} but got $observed."
    exit 1;
  fi
done
for i in $(seq 0 $lastIndex); do
  j=$((i + 1));
  #observed="$(grep albertii Efaecium_batch_example_summary.txt | cut -f $j -d $'\t')"
  observed="$(grep albertii Efaecium_batch_example_summary.txt)"
  echo "${expected3[$i]} <=> $observed";
  if [ "${expected3[$i]}" != "$observed" ]; then
    echo "ERROR: test output was incorrect in field $i.  I expected ${expected3[$i]} but got $observed."
    exit 1;
  fi
done
echo " [OK] All fields matched!"
