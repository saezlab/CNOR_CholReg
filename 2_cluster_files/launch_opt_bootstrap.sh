for i in {1..100}
do 
echo "Welcome $i times"
   bsub -oo "3_results/results_Hela_$i.o" -eo "3_results/results_Hela_$i.e" ./run_peter_pert_bootstrap_Hela.sh

   bsub -oo "3_results/results_Hek_$i.o" -eo "3_results/results_Hek_$i.e" ./run_peter_pert_bootstrap_Hek.sh

   bsub -oo "3_results/results_HepG2_$i.o" -eo "3_results/results_HepG2_$i.e" ./run_peter_pert_bootstrap_HepG2.sh

   bsub -oo "3_results/results_Huh7_$i.o" -eo "3_results/results_Huh7_$i.e" ./run_peter_pert_bootstrap_Huh7.sh

done

