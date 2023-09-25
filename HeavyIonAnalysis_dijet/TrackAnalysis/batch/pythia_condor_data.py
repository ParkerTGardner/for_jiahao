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
for i in range(48):
   # with open("all_data_list/list_cor_{:03d}".format(i)) as fs:
   #    line_count = 0
   #    for line in fs:
   #       line_count += 1

   # line_count_frac = line_count  # for all pthatpythia
   # file_count = math.ceil(line_count_frac)

   # for j in range(1, int(file_count)):
      temp = '''
   Arguments  = {:03d}
   Queue
      '''.format(j)
      command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
