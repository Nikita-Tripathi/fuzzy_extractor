# IMPORTS
from math import ceil
import numpy as np
from secrets import randbits
from subprocess import check_output
from galois_field import GFpn
import galois
from PIL import Image
from multiprocessing import Pool, Process
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys, time, random
import mx_par



# RUNS ON PYTHON 3.8.9
#   (on my system its in /usr/bin/python3, with /usr/bin/pip3 package manager)

class FuzzyExtractor:

    # Specify subsample size, number of lockers. and LPN encryption error rate
    # Default parameters are k=43, l=10^6, err=0.12, etc.
    def __init__(self, k=43, l=1000000, err=0.11, ecc_len=1224, lbd=128, xi = 128, t=12):
        self.k = k
        self.l = l
        self.error_rate = err
        self.nu = ecc_len
        self.lbd = lbd
        self.xi = xi
        # (1-gamma)/2 = err
        self.gamma = 1 - (err * 2)
        self.dec_check = ecc_len * (0.5 - self.gamma/3)

        self.t = t
        self.ecc_msg_len = t + lbd * 2 + xi
        # Below are the accumulators for the ciphertext and the positions list
        self.ctexts = [[] for _ in range(self.l)]
        # self.positions = [[] for _ in range(self.l)]
        self.lpn_matrices = [ np.array([self.bitarr(k) for a in range(ecc_len)]) for _ in range(self.l)]
        # At l = 10^6, lpn_matrix generation will take several hours (serial code - could be less in parallel, won't be less than an hour though) 
        #   at l=10, it took 0.345 sec
        #   at l=100, 0.627 sec
        #   at l=1000, 5.369 - 5.398 sec
        #   Could be more efficient to pre-compute these on university servers in parallel (5-10 sets of 10^6 matrices) and pick them randomly
        # Irreducible polynomial for GF(2^128)
        self.irreducible_poly = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 1, 1, 1
        ]
        self.gf = GFpn(2,self.irreducible_poly)


    def bitarr(self, i):
        r = randbits(i)
        return [int(d) for d in bin(r)[2:].zfill(i)]
    
    def LPN_enc(self, A, key, msg):
        # multiply LPN matrix by the key
        d = np.matmul(A, key) % 2
        # Call LDPC code (C) and do:
        #   * Generate an error correcting code (ECC) of the msg
        #   * Add self.error_rate errors (using transmit function)
        #   * Return the "corrupted" code of msg + d (addition mod 2)
        noisy_msg = check_output(["bash", "ldpc_enc.bash", f"-s {randbits(32)}", f"-m {msg}", f"-e {self.error_rate}"])
        noisy_msg= noisy_msg.decode('ASCII').split('\n')[-2]
        # Encode noisy ecc 
        m = np.array([int(b) for b in noisy_msg])
        # Step (v): Compute ciphertext
        ctxt = m ^ d
        return ctxt

    def LPN_dec(self, A, key, ctxt):
        # multiply LPN matrix by the key
        d = np.matmul(A, key) % 2
        # Add (mod 2) d and ctxt (assuming ctxt is a numpy array)
        temp = d ^ ctxt
        # encode temp into a bitstring
        tmp = ''
        for i in temp:
            tmp += str(i)
        # Call LDPC code (C) and decode temp
        decoded = check_output(["bash", "ldpc_dec.bash", f"-c {tmp}", f"-e {self.error_rate}"])
        if (decoded.decode('ASCII').split('\n')[-3].split()[3]) == '0':
            return np.array([]) # Invalid decryption
        decoded= decoded.decode('ASCII').split('\n')[-2]
        g_i = np.array([int(b) for b in decoded])
        # Check if hamming weigth of g_i is more than some value (depends on self.error_rate)
        # if sum(temp ^ g_i) > self.dec_check: # TODO!!!! FIGURE OUT HOW SUM WILL WORK WITH LIST OF ARRAYS
        #     return np.array([]) # This is our "error" output
        return g_i
    
    def m(self, fm, fx, fL):
        t1 = time.time()
        # accumulator for sum
        acc = self.gf.elm([0])
        # accumulator for x^i
        x = fx ** 0
        for i in range(1,fL):
            x = x * fx
            acc = acc + (x * fm[i])
        print(f"Serial m took {time.time() - t1} seconds")
        return acc
    


    # w is the input string (iris encoding), n is the iris mask (used for sampling unmasked bits)
    def gen(self, w, n):
        # step 1: generate R and R_1
        R = bin(randbits(self.xi))[2:].zfill(self.xi)
        R_1 = bin(randbits(self.lbd * 2))[2:].zfill(self.lbd * 2)

        to_enc = (R + R_1).zfill(self.ecc_msg_len)

        # Generate l sets of unmasked positions 
        self.positions = [ random.SystemRandom().sample(np.flatnonzero(n).tolist(), k=self.k) for _ in range(self.l) ]
        # print(self.positions)
        # step 2: start a loop
        for i in range(self.l):
            # Get a sample of w at positions
            sample_i = np.array([w[pos] for pos in self.positions[i]])

            # LPN encryption of (R|R_1)
            p_i = self.LPN_enc(self.lpn_matrices[i], sample_i, to_enc)

            self.ctexts[i] = p_i

        # step 3:
        self.L = ceil(self.l * self.nu / self.lbd) + 4
        # step 4: Encode (concat. of all ciphertexts) to a vector m of length L-4 in GF(2^self.lbd)
        # NOTE: each entry of vector m will have length self.lbd (128)
        bigctxt = ''
        for c in self.ctexts:
            for i in c:
                bigctxt += str(i)
        m = [self.gf.elm([int(j) for j in bigctxt[i:i+self.lbd]]) for i in range(0, len(bigctxt), self.lbd)]
        # print(bigctxt)
        # print(m)
        # step 5: Parse R_1 into x and y (both have length self.lbd)
        # print(len(m), self.L, len(R_1[:self.lbd]), len(R_1[self.lbd:]))
        x = self.gf.elm([int(j) for j in R_1[:self.lbd]])
        y = self.gf.elm([int(j) for j in R_1[self.lbd:]])

        # step 6: Compute T = x^L + x^2 * m(x) + x*y
        # https://asecuritysite.com/principles/gf <--- use this for polynomial mult in GF
        # mx_ = mx_par.mx_parallel(m, x, self.L - 4)
        # print(mx_) #FIXME this does not work at all! Workers keep calling the main function over and over
        mx = self.m(m, x, self.L-4) # NOTE: this operation takes 3 seconds for l = 10, L-4 = 96
        print(mx)
        self.T = x ** self.L + ((x ** 2) * mx) + x * y
        # step 7: Output key R, self.ctexts, and self.T

        return R
        
    # w_ is W' in the paper, ciphertext and T can be found in self.ctexts and self.T respectively
    def rep(self, w_):
        for i in range(self.l):
            # Get a sample of w' at positions
            sample_i = np.array([w_[pos] for pos in self.positions[i]])

            dec = self.LPN_dec(self.lpn_matrices[i], sample_i, self.ctexts[i])
            # STEP iv
            if not (len(dec) == 0 or dec[:self.t].any()): # i.e., if dec is not None
                R, R_1 = dec[self.t:self.t + self.xi], dec[self.t + self.xi:]
                # print(R,R_1)
                bigctxt = ''
                for c in self.ctexts:
                    for i in c:
                        bigctxt += str(i)
                m = [self.gf.elm([int(j) for j in bigctxt[i:i+self.lbd]]) for i in range(0, len(bigctxt), self.lbd)]
                x = self.gf.elm([j for j in R_1[:self.lbd]])
                y = self.gf.elm([j for j in R_1[self.lbd:]])

                mx = self.m(m, x, self.L-4) # NOTE: this operation takes 3 seconds for l = 10, L-4 = 96

                T_rep = x ** self.L + ((x ** 2) * mx) + x * y
                print(self.T)
                print(T_rep)
                if (T_rep - self.T) == self.gf.elm([]): #FIXME: some errors in ==
                    print("Check passed")
                    return ''.join([str(b) for b in R])
            else:
                print("Bad decryption :", dec[:15])
        # Could use secrets.compare_digest(a, b)
        # to compare T and T' ??
        print("Checks failed")
        return None
    
    def test(self):
        a = self.gf.elm([
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 1, 1, 1
        ])
        b = self.gf.elm([0, 1, 0])

        print(a+b)
        print(a * b)
        print(a ** 5)


# IMG_OPENER 
# INPUTS:
# path: String - Path to an image to be opened
# mask: Boolean - Optional flag to indicate whether the image is a mask
# OUTPUT: Flattened (1D) NumPy array of bits.
# Depending on `mask`, this is either just one array of length 512*64 or 6 of them
def img_opener(path, mask=False):
    image = Image.open(path)
    data = np.asarray(image) % 2 # map everything to 0 or 1 since it's B&W values (0 or 255)
    if mask:
        return data.flatten()
    else:
        return [data[i:i+64].flatten() for i in range(0, 384, 64)]


def main():
    mask1 = "./test_msk/04560d632_mano.bmp"
    code1 = "./test_code/04560d632_code.bmp"
    mask2 = "./test_msk/04560d634_mano.bmp"
    code2 = "./test_code/04560d634_code.bmp"

    m1 = img_opener(mask1, mask=True)
    c1 = [ m1 & c for c in img_opener(code1) ] # XOR all 6 codes (one per Gabor filter pair) with mask here
    m2 = img_opener(mask2, mask=True)
    c2 = [ m2 & c for c in img_opener(code2) ] # XOR all 6 codes (one per Gabor filter pair) with mask here

    # print(c1)

    # Matrix generation works (but takes too long)
    # LDPC parts (bash scripts) work (may need to edit to support multithreading, but unlikely)
    # LPN Enc & Dec work (small things to add in Dec)
    # Tag calculation works (as efficient as I could get it for now)
    # Sampling - works
    t1 = time.time()
    fe = FuzzyExtractor(l=7000)
    t2 = time.time()
    print(f"Initialized (generated lpn arrays & GF(2^128)) in {t2 - t1} seconds")

    # HACK ACCORDING TO FULLERS PAPER (SECTION 4), TRANSFORM #5 HAS THE BEST RATE FOR IMAGES OF SAME IRIS
    a = fe.gen(c1[5], m1)
    t3 = time.time()
    print(f"Ran GEN in {t3 - t2} seconds")
    # # print(fe.T)
    b = fe.rep(c2[5])
    print(f"Ran REP in {time.time() - t3} seconds")

    print(a)
    print(b)


    print("no problems so far")

main()

