import ctypes, numpy as np
lib = ctypes.CDLL("./libbench_hmat.so")
lib.bench_hmat.restype = ctypes.c_int
n = 2000; dim = 1
pts = (ctypes.c_double * (n*dim))(*((np.linspace(0,1,n, endpoint=False)).tolist()))
b = (ctypes.c_double*n)(*(np.cos(2*np.pi*np.arange(n)*13.0/n).tolist()))
x = (ctypes.c_double*n)()
t_ass = ctypes.c_double(); t_fact = ctypes.c_double(); t_solve = ctypes.c_double()
logdet = ctypes.c_double(); comp = ctypes.c_size_t(); uncomp = ctypes.c_size_t(); full = ctypes.c_size_t()
rc = lib.bench_hmat(n, dim, pts, ctypes.c_double(0.01), ctypes.c_double(1e-8), ctypes.c_double(1e-6), 250, 1, 1, 3, b, x,
    ctypes.byref(t_ass), ctypes.byref(t_fact), ctypes.byref(t_solve), ctypes.byref(logdet),
    ctypes.byref(comp), ctypes.byref(uncomp), ctypes.byref(full))
print("rc", rc, "t_ass", t_ass.value, "t_fact", t_fact.value, "t_solve", t_solve.value, "logdet", logdet.value, comp.value, uncomp.value)
