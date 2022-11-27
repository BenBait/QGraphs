echo ""
echo "arg = 1 ==> check that dag(Laplacian) == Laplacian"
echo "arg = 2 ==> get success probability"
echo ""

# check for 1 arg to indicate build
# don't forget spaces around inside of [...]
if [ $# -ne 1 ];
then
    echo "Pass one argument"
    exit
fi

rm a.out
[ $1 -eq 1 ] && set -x \
             && gcc get_sparse_laplacian.c csparse/csparse.c -lm \
             && set +x 
 
[ $1 -eq 2 ] && set -x \
             && gcc -O3 get_success_prob.c csparse/csparse.c -lm \
             && set +x 

echo ""
echo "Run ./a.out with args to get desired check"
