import time
time.sleep(2)        # exceeds the test's 0.5s timeout -> timed_out / tier=slow
assert True
