from multiprocessing import Process, cpu_count
import subprocess

command = ["run.mac", "run1.mac"]

#multiprocessing.cpu_count()
def MC(i):
	subprocess.call("./scifi " + command[i], shell =True)

# Create node-local processes
shared_processes = []
#for i in range(cpu_count()):
for i in range(len(command)):
    p = Process(target= MC, args = (i,))
    shared_processes.append(p)

# Start processes
for sp in shared_processes:
    sp.start()

# Wait for all processes to finish
for sp in shared_processes:
    sp.join()
