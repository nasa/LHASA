import ee
import random
from time import sleep

def connect(max_tries=15):
    limit = 0
    while limit < max_tries:
        try:
            ee.Initialize()
            print("\nEarth Engine initialized. \n")
            return
        except:
            limit += 1
            if limit<max_tries:
                print("\nEarth Engine failed to connect " +str(limit)+ " time(s). " \
                "After "+str(max_tries)+" times this process will quit.")
                retry_in = random.uniform(10.2, 120.1)
                print('Retrying in', retry_in, "seconds.")
                sleep(retry_in) # mimic sporadic retry behavior
                
        
    print("Earth Engine connection failed.")

if __name__ == "__main__":
    connect()