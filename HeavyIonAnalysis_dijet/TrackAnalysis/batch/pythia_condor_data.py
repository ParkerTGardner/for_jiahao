import math
import subprocess
name = "pythia_batch_data_output.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = pythia_batch_data_output/pythia_sub.sh
Arguments  = 000
Log        = pythia_batch_data_output/log/submit_v0.$(Process).log
Output     = pythia_batch_data_output/out/submit_v0.$(Process).out
Error      = pythia_batch_data_output/err/submit_v0.$(Process).err
+MaxRuntime =60000
Queue
'''
for i in range(1,48):
   temp = '''
Arguments  = %03d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
