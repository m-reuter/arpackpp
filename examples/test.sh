#! /bin/sh

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

complex=false

failed=0
passed=0

testEquality() {
  # Check if file is executable
  #echo $1 | awk '{print $1}'
  fn=`echo $1 | awk '{print $1;}'`
  if [ ! -x "$fn" ]; then
    echo "[${RED}FAIL${NC}] $fn does not exist or no executable"
    ((failed++))
    return 1
  fi
  
  # Check if it runs without error
  $1 | cat /dev/null
  retstate=$?
  if [ ! "$retstate" == "0" ]; then
    echo "[${RED}FAIL${NC}] $1 non zero return state: $retstate"
    ((failed++))
    return 1
  fi
  
  # Grab result
  results=`$1 2>/dev/null | grep "$3" | awk -v temp=$4 '{print $temp}'`
  #echo "results: " $results
  if (( $(grep -c . <<<"$results") > 1 )); then
    while IFS= read -r line; do 
      arr+=( "$line" ); 
    done <<< "$results"
    #echo "arr: " $arr
    result=${arr[${5:-0}]}
  else
    result=$results
  fi
  #echo "result: " $result
  
  # Grab potential error message
  errormsg=`$1 2>&1 1>/dev/null`
  
  # Fail on error message
  if [ ! -z "$errormsg" ]; then
    echo "[${RED}FAIL${NC}] $1 TEST ERROR: \n $errormsg"
    ((failed++))
    return 1
  fi

  # Check if results is as expected or not
  if [ "$result" == "$2" ]; then
    echo "[${GREEN}OK${NC}]   $1 TEST $result"
    ((passed++))
    return 0
  else
    echo "[${RED}FAIL${NC}] $1 TEST $result != $2 EXPECTED"
    ((failed++))
    return 1
  fi
}

#pre="./"

# band
pre="./band/"

if [ "$complex" == true ]; then
  # band/complex
  #crashes on my mac, or no eval output
  testEquality "${pre}bcompgre" "???" "lambda\[4\]:" "2"
  testEquality "${pre}bcompgsh" "???" "lambda\[4\]:" "2"
  testEquality "${pre}bcompreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}bcompshf" "???" "lambda\[4\]:" "2"
fi

# band/nonsym
testEquality "${pre}bnsymgre" "120979" "lambda\[4\]:" "2"
testEquality "${pre}bnsymgsc" "0.361936" "lambda\[4\]:" "2"
testEquality "${pre}bnsymgsh" "182.929" "lambda\[4\]:" "2"
testEquality "${pre}bnsymreg" "280.417" "lambda\[4\]:" "2"
testEquality "${pre}bnsymshf" "868.92" "lambda\[4\]:" "2"
testEquality "${pre}bsvd" "0.0225147" "sigma \[4\]:" "3"

# band/sym
testEquality "${pre}bsymgbkl" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}bsymgcay" "356.338" "lambda\[4\]:" "2"
testEquality "${pre}bsymgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}bsymgshf" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}bsymreg" "948.395" "lambda\[4\]:" "2"
testEquality "${pre}bsymshf" "76.8333" "lambda\[4\]:" "2"




# dense/complex
pre="./dense/"
if [ "$complex" == true ]; then
  # crashes on mac
  testEquality "${pre}dcompgre" "???" "lambda\[4\]:" "2"
  testEquality "${pre}dcompgsh" "???" "lambda\[4\]:" "2"
  testEquality "${pre}dcompreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}dcompshf" "???" "lambda\[4\]:" "2"
fi 

# dense/nonsym
testEquality "${pre}dnsymgre" "120979" "lambda\[4\]:" "2"
testEquality "${pre}dnsymgsc" "0.361936" "lambda\[4\]:" "2"
testEquality "${pre}dnsymgsh" "182.929" "lambda\[4\]:" "2"
testEquality "${pre}dnsymreg" "280.417" "lambda\[4\]:" "2"
testEquality "${pre}dnsymshf" "868.92" "lambda\[4\]:" "2"
testEquality "${pre}dsvd" "0.557234" "sigma \[4\]:" "3"
# matrix.dat in same dir :
testEquality "${pre}dsvd2" "1.36095" "sigma \[4\]:" "3"

# dense/sym
testEquality "${pre}dsymgbkl" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}dsymgcay" "356.338" "lambda\[4\]:" "2"
testEquality "${pre}dsymgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}dsymgshf" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}dsymreg" "948.395" "lambda\[4\]:" "2"
testEquality "${pre}dsymshf" "76.8333" "lambda\[4\]:" "2"




# product/complex
pre="./product/"
# crashes on mac
if [ "$complex" == true ]; then
  testEquality "${pre}compgreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}compgshf" "???" "lambda\[4\]:" "2"
  testEquality "${pre}compreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}compshf" "???" "lambda\[4\]:" "2"
fi 

# product/nonsym
# ATTENTION, the float part is failing 
testEquality "${pre}nsymgreg" "20163.1" "lambda\[4\]:" "2"
testEquality "${pre}nsymgsci" "0.5" "lambda\[4\]:" "2"
testEquality "${pre}nsymgscr" "0.5" "lambda\[4\]:" "2"
testEquality "${pre}nsymgshf" "30.4882" "lambda\[4\]:" "2"
testEquality "${pre}nsymreg" "919.781" "lambda\[4\]:" "2"
testEquality "${pre}nsymshf" "0.912926" "lambda\[4\]:" "2"
testEquality "${pre}svd" "0.0410123" "sigma\[4\]:" "2"

# product/sym
testEquality "${pre}symgbklg" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}symgcayl" "356.338" "lambda\[4\]:" "2"
testEquality "${pre}symgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}symgshft" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}symreg" "76.8333" "lambda\[4\]:" "2"
testEquality "${pre}symshft" "???" "lambda\[4\]:" "2"
testEquality "${pre}symsimp" "157.71" "lambda\[4\]:" "2"


# reverse
# crashes on mac
pre="./reverse/"
if [ "$complex" == true ]; then
  testEquality "${pre}rcompgre" "???" "lambda\[4\]:" "2"
  testEquality "${pre}rcompgsh" "???" "lambda\[4\]:" "2"
  testEquality "${pre}rcompreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}rcompshf" "???" "lambda\[4\]:" "2"
fi 

testEquality "${pre}rnsymgre" "???" "lambda\[4\]:" "2"
testEquality "${pre}rnsymgsc" "???" "lambda\[4\]:" "2"
testEquality "${pre}rnsymgsh" "???" "lambda\[4\]:" "2"
testEquality "${pre}rnsymreg" "???" "lambda\[4\]:" "2"
testEquality "${pre}rnsymshf" "???" "lambda\[4\]:" "2"
testEquality "${pre}rsvd" "???" "sigma\[4\]:" "2"

testEquality "${pre}rsymgbkl" "???" "lambda\[4\]:" "2"
testEquality "${pre}rsymgcay" "???" "lambda\[4\]:" "2"
testEquality "${pre}rsymgreg" "???" "lambda\[4\]:" "2"
testEquality "${pre}rsymgshf" "???" "lambda\[4\]:" "2"
testEquality "${pre}rsymreg" "???" "lambda\[4\]:" "2"
testEquality "${pre}rsymshf" "???" "lambda\[4\]:" "2"


# SUPERLU Stuff

# areig
pre="./areig/"
if [ "$complex" == true ]; then
  testEquality "${pre}acompgre" "???" "lambda\[4\]:" "2"
  testEquality "${pre}acompgsh" "???" "lambda\[4\]:" "2"
  testEquality "${pre}acompreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}acompshf" "???" "lambda\[4\]:" "2"
  # from non sym dir, but also complex
  testEquality "${pre}ansymgsc" "???" "lambda\[4\]:" "2"
fi 

testEquality "${pre}ansymgre" "120979" "lambda\[4\]:" "2"
testEquality "${pre}ansymgsh" "182.929" "lambda\[4\]:" "2"
testEquality "${pre}ansymreg" "919.781" "lambda\[4\]:" "2"
testEquality "${pre}ansymshf" "-0.307351" "lambda\[4\]:" "2"
testEquality "${pre}simple" "402.193" "lambda\[4\]:" "2"

testEquality "${pre}asymgbkl" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}asymgcay" "356.338" "lambda\[4\]:" "2"
testEquality "${pre}asymgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}asymgshf" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}asymreg" "76.8333" "lambda\[4\]:" "2"
testEquality "${pre}asymshf" "157.71" "lambda\[4\]:" "2"



# superlu
pre="./superlu/"
if [ "$complex" == true ]; then
  testEquality "${pre}lcompgre" "???" "lambda\[4\]:" "2"
  testEquality "${pre}lcompgsh" "???" "lambda\[4\]:" "2"
  testEquality "${pre}lcompreg" "???" "lambda\[4\]:" "2"
  testEquality "${pre}lcompshf" "???" "lambda\[4\]:" "2"
fi 

testEquality "${pre}lnsymgre" "120979" "lambda\[4\]:" "2"
# strange, this is complex, but seems to work:
testEquality "${pre}lnsymgsc" "0.5" "lambda\[4\]:" "2"
testEquality "${pre}lnsymgsh" "182.929" "lambda\[4\]:" "2"
testEquality "${pre}lnsymreg" "919.781" "lambda\[4\]:" "2"
testEquality "${pre}lnsymshf" "-0.307351" "lambda\[4\]:" "2"
testEquality "${pre}lsvd" "1.41684" "smallest singular" "4" "0"
testEquality "${pre}lsvd2" "9.89757" "sigma\[5\]:" "2"

testEquality "${pre}lsymgbkl" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}lsymgcay" "356.338" "lambda\[4\]:" "2"
testEquality "${pre}lsymgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}lsymgshf" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}lsymreg" "48.2193" "lambda\[2\]:" "2"
testEquality "${pre}lsymshf" "157.71" "lambda\[4\]:" "2"


#harwell
pre="./harwell/"
if [ "$complex" == true ]; then
  testEquality "${pre}hcompgen" "???" "lambda\[4\]:" "2"
  testEquality "${pre}hcompstd" "???" "lambda\[4\]:" "2"
fi 
testEquality "${pre}hnsymgen ${pre}mhd416a.rua ${pre}mhd416b.rua" "-1.29284" "lambda\[6\]:" "2"
testEquality "${pre}hnsymstd ${pre}mhd416a.rua" "-0.042968" "lambda\[4\]:" "2"
testEquality "${pre}hsymgen ${pre}lund_a.rsa ${pre}lund_b.rsa" "2.20462e+06" "lambda\[5\]:" "2"
testEquality "${pre}hsymstd ${pre}lund_a.rsa" "2.23854e+08" "lambda\[5\]:" "2"


# umfpack
pre="./umfpack/"
testEquality "${pre}usymgbkl" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}usymgcay" "356.338" "lambda\[4\]:" "2"
testEquality "${pre}usymgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}usymgshf" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}usymreg" "76.8333" "lambda\[4\]:" "2"
testEquality "${pre}usymshf" "157.71" "lambda\[4\]:" "2"


# cholmod
pre="./cholmod/"
testEquality "${pre}csymgreg" "122323" "lambda\[4\]:" "2"
testEquality "${pre}csymgshf" "158.117" "lambda\[4\]:" "2"
testEquality "${pre}csymreg" "76.8333" "lambda\[4\]:" "2"
testEquality "${pre}csymshf" "157.71" "lambda\[4\]:" "2"


echo
echo PASSED: $passed
echo FAILED: $failed
echo



