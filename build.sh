echo ""
echo "arg = 0 ==> build matrix exponential library"
echo "arg = 1 ==> check that dag(Laplacian) == Laplacian"
echo "arg = 2 ==> get single success probability, with parameters:"
echo "n, p, w, t, gamma"
echo ""
echo "arg = 3 ==> get success probability for spectrum given by parameters:"
echo "n, p, w, t_max, gamma_max, gamma_step"
echo ""
echo "arg = 4 ==> check eigenvectors of complete graph"
echo ""

# check for 1 arg to indicate build
# don't forget spaces around inside of [...]
if [ $# -ne 1 ];
then
    echo "Pass one argument"
    exit
fi

rm a.out
[ $1 -eq 0 ] && set -x \
             && gcc -O2 -c exp/matrix_exponential.c exp/r8lib.c exp/c8lib.c \
             && set +x 

[ $1 -eq 1 ] && set -x \
             && gcc get_sparse_laplacian.c -lm \
             && set +x 
 
[ $1 -eq 2 ] && set -x \
             && gcc -O2 -c get_success_prob.c -lm \
             && gcc -O2 get_success_prob.o matrix_exponential.o r8lib.o c8lib.o -lm \
             && set +x 

[ $1 -eq 3 ] && set -x \
             && gcc -O2 -c get_success_prob_map_data.c -lm \
             && gcc -O2 get_success_prob_map_data.o matrix_exponential.o r8lib.o c8lib.o -lm \
             && set +x 

 #&& gcc -O3 get_success_prob.c csparse/csparse.c -lm \
[ $1 -eq 4 ] && set -x \
             && gcc -O3 check_complete_eigvectors.c -lm \
             && set +x 

echo ""
echo "Run ./a.out with args to get desired check"
