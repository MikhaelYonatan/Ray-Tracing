from math import sqrt
from PIL import Image

# ==========================================
# PENDEKATAN VEKTOR EUCLIDEAN
# ==========================================
def dot_euc(a, b): 
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def sub_euc(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def mul_euc(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def reflection_euclidean(D, N):
    """R = D - 2(D·N)N"""
    return sub_euc(D, mul_euc(N, 2 * dot_euc(D, N)))


# ==========================================
# PENDEKATAN ALJABAR GEOMETRI
# ==========================================
class Multivector:
    """
    Representasi Multivektor sederhana untuk 3D.
    Hanya mendukung Grade 0 (Skalar) dan Grade 1 (Vektor).
    Bivektor disimpan sebagai tupel (e12, e23, e31).
    """
    def __init__(self, s=0.0, v=(0.0, 0.0, 0.0), b=(0.0, 0.0, 0.0)):
        self.s = s  # Grade 0 (Skalar)
        self.v = v  # Grade 1 (Vektor)
        self.b = b  # Grade 2 (Bivektor: e12, e23, e31)

    @staticmethod
    def geometric_product(A, B):
        """Implementasi ab = a.b + a^b"""
        # Bagian Skalar (Inner Product)
        s_new = A.v[0]*B.v[0] + A.v[1]*B.v[1] + A.v[2]*B.v[2]
        
        # Bagian Bivektor (Outer Product)
        b_new = (
            A.v[0]*B.v[1] - A.v[1]*B.v[0],  # e12
            A.v[1]*B.v[2] - A.v[2]*B.v[1],  # e23
            A.v[2]*B.v[0] - A.v[0]*B.v[2]   # e31
        )
        return Multivector(s=s_new, b=b_new)

    @staticmethod
    def gp_mv_vec(mv, vec):
        """Perkalian Multivektor (s + b) dengan Vektor (v)"""
        s, b = mv.s, mv.b
        v = vec
        
        v_new = (
            s*v[0] + b[0]*v[1] - b[2]*v[2],
            s*v[1] - b[0]*v[0] + b[1]*v[2],
            s*v[2] - b[1]*v[1] + b[2]*v[0]
        )
        return Multivector(v=v_new)

def reflection_geometric_algebra(D_vec, N_vec):
    """
    Implementasi R = - N D N
    Menggunakan operasi 'Sandwich Product'
    """
    # Hitung ND (Hasilnya adalah Skalar + Bivektor)
    nd = Multivector.geometric_product(Multivector(v=N_vec), Multivector(v=D_vec))
    
    # Hitung (ND)N (Hasilnya kembali menjadi Vektor)
    ndn = Multivector.gp_mv_vec(nd, N_vec)

    # Hasil akhirnya adalah -NDN
    return (-ndn.v[0], -ndn.v[1], -ndn.v[2])


# ==========================================
# 3. SISTEM RENDER & INTERSEKSI
# ==========================================
def normalize(v):
    l = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return (v[0]/l, v[1]/l, v[2]/l)

sphere_center = (0,0,3)
sphere_radius = 1.0
sphere2_center = (1, 1, 1)
sphere2_radius = 0.1
light_pos = (2, 2, 0)
camera = (0,0,0)

def intersect_sphere(ray_o, ray_d, center, radius):
    oc = sub_euc(ray_o, center)
    b = 2 * dot_euc(oc, ray_d)
    c = dot_euc(oc, oc) - radius**2
    disc = b*b - 4*c
    if disc < 0: return None
    t = (-b - sqrt(disc)) / 2
    return t if t > 0 else None

def intersect_scene(ray_o, ray_d):
    """Cek interseksi dengan semua objek, kembalikan t terdekat"""
    t1 = intersect_sphere(ray_o, ray_d, sphere_center, sphere_radius)
    t2 = intersect_sphere(ray_o, ray_d, sphere2_center, sphere2_radius)
    candidates = [t for t in [t1, t2] if t is not None]
    return min(candidates) if candidates else None

def compute_color(q):
    N = normalize(sub_euc(q, sphere_center))
    L = normalize(sub_euc(light_pos, q))
    
    D = normalize(sub_euc(q, camera))

    # PERBANDINGAN
    shadow_t = intersect_scene(q, L)
    base = 0.2 if shadow_t else max(0.2, dot_euc(N, L))

    R1 = reflection_euclidean(D, N)
    R2 = reflection_geometric_algebra(D, N)

    gloss1 = max(0, dot_euc(normalize(R1), L))
    gloss2 = max(0, dot_euc(normalize(R2), L))

    col_eu = (int(255 * (base + 0.5*gloss1)), 80, 80)
    col_ga = (int(255 * (base + 0.5*gloss2)), 80, 80)
    return col_eu, col_ga

# PROSES RENDERING
W, H = 400, 300
img = Image.new("RGB", (W*2, H))
pixels = img.load()

for y in range(H):
    for x in range(W):
        xx, yy = (x - W/2) / W, (y - H/2) / H
        ray_d = normalize((xx, yy, 1))
        t = intersect_scene(camera, ray_d)
        
        if t is None:
            pixels[x, y] = pixels[x+W, y] = (30, 30, 40)
            continue

        q = (camera[0]+ray_d[0]*t, camera[1]+ray_d[1]*t, camera[2]+ray_d[2]*t)
        col_eu, col_ga = compute_color(q)
        pixels[x, y], pixels[x+W, y] = col_eu, col_ga

img.save("raytrace_final.png")
print("Berhasil merender! Cek file raytrace_final.png")