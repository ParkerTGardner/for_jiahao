import math
import subprocess
name = "pythia_batch_output.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = pythia_batch_output/pythia_sub.sh
Arguments  = 000
Log        = pythia_batch_output/log/submit_v0.$(Process).log
Output     = pythia_batch_output/out/submit_v0.$(Process).out
Error      = pythia_batch_output/err/submit_v0.$(Process).err
+MaxRuntime =10000
Queue
'''

with open('../fileLists/test.txt') as fs:
    line_count = 0
    for line in fs:
        line_count += 1

line_count_frac = line_count/20.0 #for all pthatpythia
file_count = math.ceil(line_count_frac)

for i in range(1, int(file_count)):
   temp = '''
Arguments  = %03d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
