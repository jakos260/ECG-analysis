from qtripy import QTripy
import time

if __name__ == "__main__":
    q = QTripy()
    q.begin(1041)

    print("QTripy server running...")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        q.close()
