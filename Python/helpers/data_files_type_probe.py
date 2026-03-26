import sys, os
import numpy as np

FN = sys.argv[1] if len(sys.argv) > 1 else './PPD1_ECGSIM_NEW/PPD1_ECGSIM_NEW_ventricles.heartdist'
N = 10

if not os.path.exists(FN):
    print("File not found:", FN); sys.exit(2)

b = open(FN,'rb').read()
print("file:", FN, "size:", len(b))

if len(b) >= 8:
    hdr = b[:8]
    r0_le = int.from_bytes(hdr[0:4], 'little', signed=False)
    r1_le = int.from_bytes(hdr[4:8], 'little', signed=False)
    r0_be = int.from_bytes(hdr[0:4], 'big', signed=False)
    r1_be = int.from_bytes(hdr[4:8], 'big', signed=False)
    print("header uint32 (LE):", r0_le, r1_le, "product", r0_le*r1_le)
    print("header uint32 (BE):", r0_be, r1_be, "product", r0_be*r1_be)
else:
    print("file < 8 bytes, cannot inspect header")

payload = b[8:] if len(b) >= 8 else b
cands = [
    ('>f8','be',8),
    ('<f8','le',8),
    ('>f4','be',4),
    ('<f4','le',4),
]

print("\nTry interpreations (first %d values + stats):" % N)
for dtype, label, bytes_per in cands:
    try:
        arr = np.frombuffer(payload, dtype=dtype)
        if arr.size == 0:
            print(f"{dtype}: empty")
            continue
        # if header seemed valid, try reshape to header dims
        shape_info = ""
        if len(b) >= 8:
            if label == 'le':
                r, c = r0_le, r1_le
            else:
                r, c = r0_be, r1_be
            prod = r*c
            if prod == arr.size or prod <= arr.size:
                try:
                    arr2 = arr[:prod].reshape((r,c))
                    shape_info = f" -> reshaped ({r},{c})"
                    sample = arr2.flatten()[:N]
                    stats = (float(np.nanmin(arr2)), float(np.nanmax(arr2)), float(np.nanmean(arr2)))
                except Exception:
                    arr2 = arr
                    sample = arr2[:N]
                    stats = (float(np.nanmin(arr2)), float(np.nanmax(arr2)), float(np.nanmean(arr2)))
            else:
                sample = arr[:N]; stats=(float(np.nanmin(arr)), float(np.nanmax(arr)), float(np.nanmean(arr)))
        else:
            sample = arr[:N]; stats=(float(np.nanmin(arr)), float(np.nanmax(arr)), float(np.nanmean(arr)))
        print(f"{dtype}{shape_info}: sample={sample.tolist()} min/max/mean={stats}")
    except Exception as e:
        print(dtype, "error:", e)

# Also try to load with current project loadmat if available
try:
    from reader import loadmat
    arr = loadmat(FN)
    print("\nreader.loadmat returned shape:", np.asarray(arr).shape,
          "min/max/mean:", float(np.nanmin(arr)), float(np.nanmax(arr)), float(np.nanmean(arr)))
    print("first rows:", np.asarray(arr).flatten()[:N].tolist())
except Exception as e:
    print("\nreader.loadmat not available or failed:", e)