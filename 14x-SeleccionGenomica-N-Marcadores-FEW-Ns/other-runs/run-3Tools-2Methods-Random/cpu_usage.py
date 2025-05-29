import psutil
import time

while True:
    print("CPU Usage per Core:", psutil.cpu_percent(percpu=True))
    time.sleep(1)

