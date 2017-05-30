#CVC
##---NO ENVIRONMENT
    ## Zero
          #./Release_Intel64/gmaster13 inputfile_weak.inp computation_zerovolt
          #sleep 1 
        #  ./Release_Intel64/gmaster13 inputfile_strong.inp computation_zerovolt
        #  sleep 1 
        # ./Release_Intel64/gmaster13 inputfile_verystrong.inp computation_zerovolt_low
        # sleep 1

    # Nonzero
        #./Release_Intel64/gmaster13 inputfile_weak_nonzero.inp computation_nonzero
          #sleep 1 
        #./Release_Intel64/gmaster13 inputfile_weak_nonzero_t10.inp computation_nonzero
          #sleep 1 
        #  ./Release_Intel64/gmaster13 inputfile_strong_nonzero.inp computation_nonzero
        #  sleep 1 
        #./Release_Intel64/gmaster13 inputfile_verystrong_nonzero.inp computation_nonzero_low
        #sleep 1
        
        #Delocal
            #./Release_Intel64/gmaster13 inputfile_weak_nonzero_deloc_env_7.inp computation_nonzero_deloc_7
            #sleep 1 

            #./Release_Intel64/gmaster13 inputfile_weak_nonzero_deloc_env_8.inp computation_nonzero_deloc_8
            #sleep 1 

            #./Release_Intel64/gmaster13 inputfile_weak_nonzero_deloc_env_9.inp computation_nonzero_deloc_9
            #sleep 1 

            #Without delocalized states
            #./Release_Intel64/gmaster13 inputfile_weak_nonzero_deloc.inp computation_nonzero_deloc
            #sleep 1 
    ##Without delocalized states but with environment
    #./Release_Intel64/gmaster13 inputfile_weak_nonzero_deloc_env.inp computation_nonzero_deloc
    #sleep 1 

            #./Release_Intel64/gmaster13 inputfile_weak_nonzero_15.inp computation_nonzero
            #sleep 1 
         #./Release_Intel64/gmaster13 inputfile_verystrong_comp_gated.inp computation_zerovolt_low
         #sleep 1


    ##Without delocalized states but with environment and with U_b = 1.5 V
    #/Release_Intel64/gmaster13 inputfile_weak_nonzero_deloc_env_15.inp computation_nonzero_deloc
#sleep 1 

#### Hypothesis on the medium system
#     ./Release_Intel64/gmaster13 inputfile_strong_nonzero_2_1.inp computation_nonzero
#     sleep 1 
      #./Release_Intel64/gmaster13 inputfile_strong_nonzero_1_7.inp computation_nonzero
#     sleep 1 
      #./Release_Intel64/gmaster13 inputfile_strong_nonzero_1_3.inp computation_nonzero
      #sleep 1 



#Hypothesis Small Barrier On Weak
##Smaller Barrier for Weak
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_sym_1.inp computation_zerovolt
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_sym_2.inp computation_zerovolt
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_sym_3.inp computation_zerovolt
#sleep 1


### Switching WEAK with different sigmas
  
# ./Release_Intel64/gmaster13 inputfile_weak_sigma_0.inp computation_zerovolt
# sleep 1 
 
# ./Release_Intel64/gmaster13 inputfile_weak_sigma_1.inp computation_zerovolt
# sleep 1 

# ./Release_Intel64/gmaster13 inputfile_weak_sigma_2.inp computation_zerovolt
# sleep 1 

# ./Release_Intel64/gmaster13 inputfile_weak_sigma_3.inp computation_zerovolt
# sleep 1 

# ./Release_Intel64/gmaster13 inputfile_weak_sigma_4.inp computation_zerovolt
# sleep 1 



#    ## Temperature T = 5 K 
#   ./Release_Intel64/gmaster13 inputfile_weak_T5.inp computation_zerovolt
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_strong_T5.inp computation_zerovolt
#  sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_verystrong_T5.inp computation_zerovolt_low
#   sleep 1 
#./Release_Intel64/gmaster13 inputfile_weak_nonzero_T5.inp computation_nonzero
#sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_strong_nonzero_T5.inp computation_nonzero
#   sleep 1 
#  ./Release_Intel64/gmaster13 inputfile_verystrong_nonzero_T5.inp computation_nonzero_low
#   sleep 1 

#    ### Symmetric Case
#  ./Release_Intel64/gmaster13 inputfile_verystrong_symmetric.inp computation_nonzero_low
#  sleep 1 
# ./Release_Intel64/gmaster13 inputfile_strong_symmetric.inp computation_nonzero
# sleep 1 
#  ./Release_Intel64/gmaster13 inputfile_weak_symmetric.inp computation_nonzero
#  sleep 1 

#### WEAK
# Plus Richtung
#   ### First Resonance +++
#./Release_Intel64/gmaster13 inputfile_weak_resonance.inp computation_resonances
#sleep 1 
## Hypo: Are the population equal for another coupling of the states as a vertical shift

#./Release_Intel64/gmaster13 inputfile_weak_resonance_plus.inp computation_resonances
#sleep 1 

#./Release_Intel64/gmaster13 inputfile_weak_resonance_minus.inp computation_resonances
#sleep 1 

#   ### Second Resonance +++
#./Release_Intel64/gmaster13 inputfile_weak_second_resonance.inp computation_resonances
#sleep 1 
#  ### Third Resonace +++
#./Release_Intel64/gmaster13 inputfile_weak_third_resonance.inp computation_resonances
#sleep 1 

###With Enviroment
#./Release_Intel64/gmaster13 inputfile_weak_resonance_env.inp computation_resonances
#sleep 1 

# Negative Richtung
#./Release_Intel64/gmaster13 inputfile_weak_resonance_neg.inp computation_resonances
#sleep 1 
#   ### Second Resonance +++
#./Release_Intel64/gmaster13 inputfile_weak_second_resonance_neg.inp computation_resonances
#sleep 1 
#  ### Third Resonace +++
#./Release_Intel64/gmaster13 inputfile_weak_third_resonance.inp computation_resonances
#sleep 1 

### STRONG
# ./Release_Intel64/gmaster13 inputfile_strong_resonance.inp computation_resonances
# sleep 1

#    ##----Enviroment
#./Release_Intel64/gmaster13 inputfile_weak_env.inp computation_zerovolt_env
#sleep 1
# ./Release_Intel64/gmaster13 inputfile_strong_env.inp computation_zerovolt_env
# sleep 1 
#/Release_Intel64/gmaster13 inputfile_verystrong_env.inp computation_zerovolt_low_env
#leep 1 

#./Release_Intel64/gmaster13 inputfile_weak_nonzero_env.inp computation_nonzero_env
#sleep 1 
#./Release_Intel64/gmaster13 inputfile_weak_nonzero_env_deloc.inp computation_nonzero_env_deloc
#sleep 1 

# ./Release_Intel64/gmaster13 inputfile_strong_nonzero_env.inp computation_nonzero_env
# sleep 1 
#/Release_Intel64/gmaster13 inputfile_verystrong_nonzero_env.inp computation_nonzero_low_env
#leep 1


## Environment but another wc
#/Release_Intel64/gmaster13 inputfile_weak_env_wc1.inp computation_zerovolt_env
#leep 1 
#./Release_Intel64/gmaster13 inputfile_weak_nonzero_env_wc1.inp computation_nonzero_env
#leep 1 


### Hypothesen Test
#### Quaratische Funktionen mit switch
#/Release_Intel64/gmaster13 inputfile_quadratisch.inp computation_quadratisch
#leep 1 
#/Release_Intel64/gmaster13 inputfile_quadratisch_noswitch.inp computation_quadratisch
#leep 1 

####With Enviroment
#/Release_Intel64/gmaster13 inputfile_quadratisch_shift_env.inp computation_quadratisch
#leep 1 
#/Release_Intel64/gmaster13 inputfile_quadratisch_shift_more_env.inp computation_quadratisch
#leep 1 
#/Release_Intel64/gmaster13 inputfile_quadratisch_shift_more_more_env.inip computation_quadratisch
#leep 1 

### No Enviroment
#/Release_Intel64/gmaster13 inputfile_quadratisch_shift.inp computation_quadratisch
#leep 1 
#/Release_Intel64/gmaster13 inputfile_quadratisch_shift_more.inp computation_quadratisch
#leep 1 
#/Release_Intel64/gmaster13 inputfile_quadratisch_shift_more_more.inp computation_quadratisch
#sleep 1 

#    ##---- T = 5
#   ./Release_Intel64/gmaster13 inputfile_weak_env_T5.inp computation_zerovolt_env
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_strong_env_T5.inp computation_zerovolt_env
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_verystrong_env_T5.inp computation_zerovolt_low_env
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_weak_nonzero_env_T5.inp computation_nonzero_env
#  sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_strong_nonzero_env_T5.inp computation_nonzero_env
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_verystrong_nonzero_env_T5.inp computation_nonzero_low_env
#   sleep 1

### Switching Testbed
#./Release_Intel64/gmaster13 inputfile_weak_switching_testbed.inp computation_switching_testbed
#sleep 1 

#./Release_intel64/gmaster13 inputfile_strong_switching_env_half.inp computation_switching_testbed
#sleep 1

 
 # #Switching gif
#./Release_Intel64/gmaster13 inputfile_weak_switching_gif.inp computation_switching_gif
#sleep 1 
# ./Release_Intel64/gmaster13 inputfile_strong_switching_gif.inp computation_switching_gif
#   sleep 1
# ./Release_Intel64/gmaster13 inputfile_verystrong_switching_gif.inp computation_switching_low_gif
# sleep 1

##Switching
#./Release_Intel64/gmaster13 inputfile_weak_switching.inp computation_switching
#sleep 1 

#./Release_Intel64/gmaster13 inputfile_weak_switching_env.inp computation_switching
#sleep 1

#### Hypo: How strong the tunneling elements P_1 in relation with the# 
#./Release_Intel64/gmaster13 inputfile_weak_switching_frank_plus.inp computation_switching
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_switching_frank_minus.inp computation_switching
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_switching.inp computation_switching
#sleep 1 


#./Release_Intel64/gmaster13 inputfile_weak_switching_deloc.inp computation_switching_deloc
#sleep 1 

#./Release_Intel64/gmaster13 inputfile_weak_switching_exbath.inp computation_switching
#sleep 1 


#./Release_Intel64/gmaster13 inputfile_weak_switching_env.inp computation_switching
#sleep 1

 #./Release_Intel64/gmaster13 inputfile_strong_switching.inp computation_switching
 #sleep 1 

#./Release_Intel64/gmaster13 inputfile_strong_switching_env.inp computation_switching
 #sleep 1

 #./Release_Intel64/gmaster13 inputfile_verystrong_switching.inp computation_switching_low
 #sleep 1 
#./Release_Intel64/gmaster13 inputfile_verystrong_switching_env.inp computation_switching_low
 #sleep 1

### WEAK SWITCHING zu gucken ob der Resonance Peak sich zeigt
#./Release_Intel64/gmaster13 inputfile_weak_switching_env_peak.inp computation_switching_peak
#sleep 1



#### Heatmap
#Environment here is lower eta = 0.004
#Release_Intel64/gmaster13 inputfile_weak_env_low.inp computation_switching_env_low
#sleep 1 

#    Gegentest NO SWITCH
#./Release_Intel64/gmaster13 inputfile_weak_switching_env_noswitch.inp computation_switching_noswitch
#sleep 1

#./Release_Intel64/gmaster13 inputfile_weak_switching_noswitch.inp computation_switching_noswitch
#sleep 1


### ONE ELECTRON
  ./Release_Intel64/gmaster13 inputfile_weak_one_electron.inp computation_oneelectron
  sleep 1 


#   Exaktes Ausfahren Um den Symmetrie punkt herum.  
#   ./Release_Intel64/gmaster13 inputfile_weak_switching_env_exact.inp computation_switching_exact_large
#   sleep 1
#  ./Release_Intel64/gmaster13 inputfile_weak_switching_env_exact.inp computation_switching_exact_medium
#  sleep 1
#  ./Release_Intel64/gmaster13 inputfile_weak_switching_env_exact.inp computation_switching_exact_small
#  sleep 1

#Wavefunction Switching

#./Release_Intel64/gmaster13 inputfile_weak_switching_wave.inp computation_switching_wave_large
#sleep 1 
#./Release_Intel64/gmaster13 inputfile_strong_switching_wave.inp computation_switching_wave_medium
#sleep 1 
#./Release_Intel64/gmaster13 inputfile_verystrong_switching_wave.inp computation_switching_wave_small
#sleep 1 


#    #Switching B05
#  ./Release_Intel64/gmaster13 inputfile_strong_switching_B05.inp computation_switching
#  sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_strong_switching_env_B05.inp computation_switching
#   sleep 1
#  ./Release_Intel64/gmaster13 inputfile_weak_switching_B05.inp computation_switching
#  sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_weak_switching_env_B05.inp computation_switching
#   sleep 1
#   ./Release_Intel64/gmaster13 inputfile_verystrong_switching_B05.inp computation_switching_low
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_verystrong_switching_env_B05.inp computation_switching_low
#   sleep 1

#Switching: Special
#./Release_Intel64/gmaster13 inputfile_weak_switching_T5.inp computation_switching
#sleep 1
#/Release_Intel64/gmaster13 inputfile_strong_switching_T5.inp computation_switching_special_strong
#leep 1 
#./Release_Intel64/gmaster13 inputfile_verystrong_switching_T5.inp computation_switching_low
#sleep 1 

##Switching: HypothesenTest
### Nur ein Lead
#./Release_Intel64/gmaster13 inputfile_weak_switching_T5.inp computation_switching
#sleep 1
### Nur das Enviroment

### Nächste Hpyothensen
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_1.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_2.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_3.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_env.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_1_env.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_2_env.inp computation_hypothese
#sleep 1
#./Release_Intel64/gmaster13 inputfile_weak_smaller_barrier_3_env.inp computation_hypothese
#sleep 1

#Switching - B=0.3
#   ./Release_Intel64/gmaster13 inputfile_strong_switching_03.inp computation_switching
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_weak_switching_03.inp computation_switching
#   sleep 1
#   ./Release_Intel64/gmaster13 inputfile_verystrong_switching_03.inp computation_switching_low
#   sleep 1

### Electron Affinity
#  ./Release_Intel64/gmaster13 inputfile_weak_diode_env_aff02.inp computation_diode
#  sleep 1

#  ./Release_Intel64/gmaster13 inputfile_weak_diode_env_aff03.inp computation_diode
#  sleep 1

#   ./Release_Intel64/gmaster13 inputfile_weak_diode_env_aff04.inp computation_diode
#   sleep 1

#   ./Release_Intel64/gmaster13 inputfile_weak_diode_env_aff05.inp computation_diode
#   sleep 1

#Evolution
#  ./Release_Intel64/gmaster13 inputfile_weak_resonance_switching.inp computation_evolution
#  sleep 1 
#  ./Release_Intel64/gmaster13 inputfile_weak_off_resonance_switching.inp computation_evolution
#  sleep 1 
   
#Pure 
#   ./Release_Intel64/gmaster13 inputfile_weak_resonance_pure.inp computation_evolution_pure
#   sleep 1 
#   ./Release_Intel64/gmaster13 inputfile_weak_off_resonance_pure.inp computation_evolution_pure
#   sleep 1 
#  ./Release_Intel64/gmaster13 inputfile_weak_symmetric_pure.inp computation_evolution_pure_expokit
#  sleep 1 


#TEST
# Test der neuen wavefunction-integral function
#./Release_Intel64/gmaster13 inputfile_weak_switching_test.inp computation_switching_test
#sleep 1 

#############################################################
## Frank Section i
#### Heatmap
#/Release_Intel64/gmaster13 inputfile_weak_frank.inp computation_switching_frank_heatmap
#leep 1 


#/Release_Intel64/gmaster13 inputfile_weak_frank_env.inp computation_switching_frank_heatmap
#leep 1 


#/Release_Intel64/gmaster13 inputfile_weak_frank_T5.inp computation_switching_frank_heatmap
#leep 1
#Release_Intel64/gmaster13 inputfile_weak_frank_symmetric.inp computation_switching_frank

#/Release_Intel64/gmaster13 inputfile_weak_switching.inp computation_switching_frank
#sleep 1 
#/Release_Intel64/gmaster13 inputfile_weak_switching_nodeloc.inp computation_switching_frank_nodeloc
#leep 1 

