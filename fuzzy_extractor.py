# IMPORTS
from math import ceil
import numpy as np
from secrets import randbits
from subprocess import check_output
import galois
from PIL import Image
import multiprocessing
import hashlib

import sys, time, random
import mx_par



# RUNS ON PYTHON 3.8.9 +
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
        
        # Below are the LPN matrices
        # self.lpn_matrices = [ np.array([self.bitarr(k) for a in range(ecc_len)]) for _ in range(self.l) ]
        self.lpn_matrices = mx_par.generateLPN(self.bitarr, k, ecc_len, l)
        # np.save("LPN_Arrays/test.npy", self.lpn_matrices[0])

        # Irreducible polynomial for GF(2^128)
        self.irreducible_poly = galois.primitive_poly(2, 128)

        self.hash = hashlib.sha3_512()
        self.L = ceil((self.hash.digest_size * 8) / self.lbd) + 4
        self.hash = hashlib.sha3_512().name




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
    
    # def LPN_enc_batch(self, A)
    # * precompute the noisy msg X self.l in GEN
    # * call the batch enc function with all the LPN matrices & the subsamples
    def LPN_batch_enc(self, keys, msgs):
        # Multiply LPN matrices by the LPN keys (subsamples of iris code)
        d = [np.matmul(self.lpn_matrices[i], keys[i]) % 2 for i in range(self.l)]
        # Prep the messages to encode (put each on a new line)
        # msg = '\n'.join(msgs)
        with open('src.src', 'w') as f:
            f.writelines(msgs)
        # Encode messages and parse the output appropriately
        noisy_msg = check_output(["bash", "ldpc_enc_batch.bash", f"-s {randbits(32)}", f"-m src.src", f"-e {self.error_rate}"])

        noisy_msg = noisy_msg.decode('ASCII').split('\n')[:-1]
        # Transform the str output to a binary vector
        m = [np.array([int(b) for b in nm]) for nm in noisy_msg]
        # Compute l ciphetexts
        ctxt = [m[i] ^ d[i] for i in range(self.l)]
        return ctxt

    def LPN_dec(self, A, key, ctxt, process):
        # multiply LPN matrix by the key
        d = np.matmul(A, key) % 2
        # Add (mod 2) d and ctxt (assuming ctxt is a numpy array)
        temp = d ^ ctxt
        # encode temp into a bitstring
        tmp = ''
        for i in temp:
            tmp += str(i)
        # Call LDPC code (C) and decode temp
        decoded = check_output(["bash", "ldpc_dec_batch.bash", f"-c {tmp}", f"-e {self.error_rate}", "-p", str(process)])
        # print(type(process), decoded)
        if (decoded.decode('ASCII').split()[3]) == '0':
            return np.array([]) # Invalid decryption
        # with 
        decoded= decoded.decode('ASCII').split()[-1]
        g_i = np.array([int(b) for b in decoded])
        # Check if hamming weigth of g_i is more than some value (depends on self.error_rate)
        # if sum(temp ^ g_i) > self.dec_check: # TODO!!!! FIGURE OUT HOW SUM WILL WORK WITH LIST OF ARRAYS
        #     return np.array([]) # This is our "error" output
        return g_i

    def LPN_dec_batch(self, As, keys, ctxts, process):
        # Decypts a batch of ciphertexts in one go (in one call)
        # multiply LPN matrix by the key
        # d = np.matmul(A, keys) % 2
        d = [np.matmul(As[i], keys[i]) % 2 for i in range(len(ctxts))]

        # Add (mod 2) d and ctxt (assuming ctxt is a numpy array)
        # temp = [d[i] ^ ctxts[i] for i in range(len(ctxts))]

        tmp = ''
        for i in range(len(ctxts)):
            for j in (d[i] ^ ctxts[i]):
                tmp += str(j)

        print(f'Testing LPN_dec_batch. Process id: {process}\n temp: {len(tmp)}')
        # encode temp into a bitstring
        input_file_name = f'r{process}.rec'
        with open(input_file_name, 'w') as f:
            f.write(tmp)
        # Call LDPC code (C) and decode temp
        decoded = check_output(["bash", "ldpc_dec_batch.bash", f"-c {input_file_name}", f"-e {self.error_rate}", "-p", str(process)])
        decoded = decoded.decode("ASCII").split("\n")
        details, test = [i.split() for i in decoded[1:-3]], decoded[-3]
        # print(f'Testing LPN_dec_batch. Process id: {process}\n bash script output:',details)
        if (test.split()[3]) == '0':
            print(f'Testing LPN_dec_batch. Process id: {process}. All invalid...')
            return np.array([]) # Invalid decryption
        
        print(f'Testing LPN_dec_batch. Process id: {process}\nFound valid decoding!')
        valid_msg_index = -1
        for i in details:
            # Layout of i: [index_of_msg, num_of_iterations, valid/invalid (0/1), change%]; 
            if i[2] == '1':
                valid_msg_index = int(i[0])
                break
        

        output_file_name = f'e{process}.ext'
        with open(output_file_name, 'r') as f:
            out = f.readlines()

        decoded = out[valid_msg_index].strip()
        
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
    

    def mac(self, key, ciphertexts):
        t = time.time()
        h = hashlib.new(self.hash)
        bigctxt = ''
        for c in ciphertexts:
            for i in c:
                bigctxt += str(i)
        # print(f"R_1: {key}\nCiphertext: {bigctxt}")
        # Generate a digest of entire ciphertexts
        bigctxt = bigctxt.encode()
        h.update(bigctxt)
        bigctxt = h.hexdigest()
        bigctxt = bin(int(bigctxt, base=16))[2:].zfill(512)

        # Encode digest into vector m
        m = [galois.Poly.Int(int(bigctxt[i:i+self.lbd], base=2)) for _ in range(0, len(bigctxt), self.lbd)]

        # split key into x and y
        x = galois.Poly.Int(int(key[:self.lbd], base=2))
        y = galois.Poly.Int(int(key[self.lbd:], base=2))

        mx = mx_par.mx_serial(m, x, self.L-4, self.irreducible_poly) 
        # mx = mx_par.mx_parallel(m, x, self.L-4, self.irreducible_poly) 

        T_rep = (pow(x, self.L, self.irreducible_poly) + (pow(x, 2, self.irreducible_poly) * mx) + (x * y)) % self.irreducible_poly
        
        print(f'Calculated MAC in {time.time() - t} seconds')

        return T_rep
        


    # w is the input string (iris encoding), n is the iris mask (used for sampling unmasked bits)
    def gen(self, w, n):
        # step 1: generate R and R_1
        R = bin(randbits(self.xi))[2:].zfill(self.xi)
        R_1 = bin(randbits(self.lbd * 2))[2:].zfill(self.lbd * 2)

        to_enc = (R + R_1).zfill(self.ecc_msg_len)

        # Generate l sets of unmasked positions 
        self.positions = [ random.SystemRandom().sample(np.flatnonzero(n).tolist(), k=self.k) for _ in range(self.l) ]

        samples = []
        # step 2: start a loop
        for i in range(self.l):
            # Get a sample of w at positions
            sample_i = np.array([w[pos] for pos in self.positions[i]])
            samples.append(sample_i)


        self.ctexts = self.LPN_batch_enc(samples, [to_enc for _ in range(self.l)])

        self.T = self.mac(R_1, self.ctexts) 

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
                R = ''
                for c in dec[self.t:self.t + self.xi]:
                    R += str(c)
                
                R_1 = ''
                for c in dec[self.t + self.xi:]:
                    R_1 += str(c)
                print("R_1:", R_1)

                T_rep = self.mac(R_1, self.ctexts)

                print(self.T)
                print(T_rep)

                if T_rep  == self.T:
                    print("Check passed")
                    return ''.join([str(b) for b in R])
            # else:
            #     print("Bad decryption :", dec[:15])

        print("Checks failed")
        return None

    def rep_parallel(self, w, num_processes=1):
        finished = multiprocessing.Array('b', False)
        split = np.array_split(range(self.l), num_processes)
        finished = multiprocessing.Manager().list([None for _ in range(num_processes)])
        processes = []
        for x in range(num_processes):
            p = multiprocessing.Process(
                target=self.rep_process, args=(w, split[x], finished, x)
            )
            processes.append(p)
            p.start()
        for p1 in processes:
            p1.join()
        if any(finished):
            print("Rep succeeded")
            return next(item for item in finished if item is not None)
        print("Rep failed")
        return None
    
    def rep_process(self, w_, indices, finished, process_id):
        counter = 0 # Track how many lockers we've checked
        samples = []
        matrices = []
        ctxts = []
        t1 = time.time()
        for i in indices:
            sample_i = np.array([w_[pos] for pos in self.positions[i]])
            samples.append(sample_i)
            matrices.append(self.lpn_matrices[i])
            ctxts.append(self.ctexts[i])
        print(f"Rep process {process_id}: Took {time.time() - t1} seconds to gather samples, matrices, and ctexts")
        dec = self.LPN_dec_batch(matrices, samples, ctxts, process_id)

        # print(dec)
        # dec=[]
        # TODO finish this AND test....
        # STEP iv
        if not (len(dec) == 0 or dec[:self.t].any()): # i.e., if dec is not None
            R = ''
            for c in dec[self.t:self.t + self.xi]:
                R += str(c)
            
            R_1 = ''
            for c in dec[self.t + self.xi:]:
                R_1 += str(c)

            
            print(R, R_1)

            T_rep = self.mac(R_1, self.ctexts)

            print(self.T)
            print(T_rep)

            if T_rep == self.T:
                print("Check passed")
                finished[process_id] = R
                return

            # counter += 1
            # if counter == 1000:
            #     if (not any(finished)):
            #         counter = 0
            #         print(f"Counter: 1000")
            #     else:
            #         return 
        return



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


    t1 = time.time()
    fe = FuzzyExtractor(l=1000)
    t2 = time.time()
    print(f"Initialized (generated lpn arrays & GF(2^128)) in {t2 - t1} seconds")

    # HACK ACCORDING TO FULLERS PAPER (SECTION 4), TRANSFORM #5 HAS THE BEST RATE FOR IMAGES OF SAME IRIS
    a = fe.gen(c1[5], m1)
    t3 = time.time()
    print(f"Ran GEN in {t3 - t2} seconds") # For l = 10000 = 10^4 typically takes 370 seconds
    # # print(fe.T)
    # b = fe.rep(c2[5]) # For l = 10000 = 10^4 typically takes 370 seconds
    # t4 = time.time()
    # print(f"Ran REP in {t4 - t3} seconds")
    # c = fe.rep_parallel(c2[5], num_processes=4)
    c = fe.rep_parallel(c2[5], num_processes=multiprocessing.cpu_count())
    # c = fe.rep_parallel(c2[5], num_processes=6)
    print(f"Ran REP parallel in {time.time() - t3} seconds")

    print(a)
    # print(b)
    print(c)


    print("no problems so far")


if __name__ == '__main__':

    main()

