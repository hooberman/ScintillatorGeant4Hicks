cd ../build
make
cd -
rm -f output/sipm_hits.txt


#for i in {1..24}; do

#sleep 3600

#

#for i in {1..2}; do

#for i in {49..96}; do

for i in {97..120}; do
    f="output/sipm_hits_Trial1_240kevents_10000events_run${i}.log"
    t="output/sipm_hits_Trial1_240kevents_10000events_run${i}.txt"
    
    echo "../build/sim_bc408_1 -n 10000 -l \"Trial1_240kevents\" -r ${i} --novis > ${f} 2>&1 &"

    rm -f "$f"
    rm -f "$t"

    ../build/sim_bc408_1 -n 10000 -l "Trial1_240kevents" -r "$i" --novis > "$f" 2>&1 &
    sleep 30
done

sleep 28800



for i in {121..144}; do
    f="output/sipm_hits_Trial1_240kevents_10000events_run${i}.log"
    t="output/sipm_hits_Trial1_240kevents_10000events_run${i}.txt"
    
    echo "../build/sim_bc408_1 -n 10000 -l \"Trial1_240kevents\" -r ${i} --novis > ${f} 2>&1 &"

    rm -f "$f"
    rm -f "$t"

    ../build/sim_bc408_1 -n 10000 -l "Trial1_240kevents" -r "$i" --novis > "$f" 2>&1 &
    sleep 30
done

sleep 28800

for i in {145..168}; do
    f="output/sipm_hits_Trial1_240kevents_10000events_run${i}.log"
    t="output/sipm_hits_Trial1_240kevents_10000events_run${i}.txt"
    
    echo "../build/sim_bc408_1 -n 10000 -l \"Trial1_240kevents\" -r ${i} --novis > ${f} 2>&1 &"

    rm -f "$f"
    rm -f "$t"

    ../build/sim_bc408_1 -n 10000 -l "Trial1_240kevents" -r "$i" --novis > "$f" 2>&1 &
    sleep 30
done

sleep 28800



for i in {169..192}; do
    f="output/sipm_hits_Trial1_240kevents_10000events_run${i}.log"
    t="output/sipm_hits_Trial1_240kevents_10000events_run${i}.txt"
    
    echo "../build/sim_bc408_1 -n 10000 -l \"Trial1_240kevents\" -r ${i} --novis > ${f} 2>&1 &"

    rm -f "$f"
    rm -f "$t"

    ../build/sim_bc408_1 -n 10000 -l "Trial1_240kevents" -r "$i" --novis > "$f" 2>&1 &
    sleep 30
done

sleep 28800

for i in {193..216}; do
    f="output/sipm_hits_Trial1_240kevents_10000events_run${i}.log"
    t="output/sipm_hits_Trial1_240kevents_10000events_run${i}.txt"
    
    echo "../build/sim_bc408_1 -n 10000 -l \"Trial1_240kevents\" -r ${i} --novis > ${f} 2>&1 &"

    rm -f "$f"
    rm -f "$t"

    ../build/sim_bc408_1 -n 10000 -l "Trial1_240kevents" -r "$i" --novis > "$f" 2>&1 &
    sleep 30
done

sleep 28800



for i in {217..240}; do
    f="output/sipm_hits_Trial1_240kevents_10000events_run${i}.log"
    t="output/sipm_hits_Trial1_240kevents_10000events_run${i}.txt"
    
    echo "../build/sim_bc408_1 -n 10000 -l \"Trial1_240kevents\" -r ${i} --novis > ${f} 2>&1 &"

    rm -f "$f"
    rm -f "$t"

    ../build/sim_bc408_1 -n 10000 -l "Trial1_240kevents" -r "$i" --novis > "$f" 2>&1 &
    sleep 30
done

sleep 28800









# 2:20 (140 seconds) for 100 events -> 14000 seconds ~ 4 hours for 10,000 events


#../build/sim_bc408_1 -n 1 -l "test2Rounded"
#../build/sim_bc408_1 -n 10 -l "Theta30_Phi90"
#../build/sim_bc408_1 -n 100 > output/log100.txt 2>&1 & ### 3mins
#../build/sim_bc408_1 -n 10000 > output/log10000.txt 2>&1 & ### 3mins
#../build/sim_bc408_1 -n 1000 > log100.txt 2>&1 &
#../build/sim_bc408_1 -n 1000 > log1000.txt 2>&1 &
