from sage.crypto.boolean_function import BooleanFunction

#	dim_redund = ceil(log(1/(corr_capa^2),2))
def Encoder(dim_mess,dim_redund):
	if (dim_redund > dim_mess):
		return None
	R_base = random_matrix(GF(2), dim_redund, dim_mess)
	while not (R_base.rank() == dim_redund):
		R_base = random_matrix(GF(2), dim_redund, dim_mess)
	I = identity_matrix(GF(2), dim_mess)
	C = identity_matrix(GF(2), dim_redund)
	G = None
	for mv in I.columns():
		mvv = matrix(mv).transpose()
		for chkcomb in C.column_space():
			chkv = matrix(chkcomb * R_base).transpose()
			chk = chkv + mvv
			if G == None:
				G = chk
			else:
				G = block_matrix([G,chk],nrows=1)
	return G

# seems to work, but expensive
def ListDecoder(codew, dim_mess, dim_redund):
	def thre(x):
		if (x > 2^(dim_redund-1)):
			return 1
		else:
			return 0

	I = identity_matrix(GF(2), dim_mess)
	C = identity_matrix(GF(2), dim_redund)
	List_res = []
	guessc = 0
	for guess in C.column_space(): # all possibilities for the parity-checks (blind)
		guessm = matrix(guess).transpose()
#		print guessm
		counters = [0]*dim_mess
		for b in range(dim_mess): # all message bits
			cvc = 0
#			print "====*===="
			for cv in C.column_space(): # all parity_checks for this bit (teh quadratic sub-obt)
				cvv = matrix(cv)
				curb = codew.column(b*2^dim_redund + cvc)	# the codeword bit noisy measurement
#				print curb
#				print "====**===="
				curh = cvv*guessm 								# the assumed parity check
#				print curh
#				print "====**===="
				cure = curb + curh[0]							# the noisy deduction
#				print cure
#				print "====**===="
				counters[b] += int(cure[0])
				cvc += 1
#		print counters
#		print map(thre, counters)
		print "====***===="
		List_res.append(matrix(GF(2), map(thre, counters)))
	return List_res

def FWHT_ListDecoder(codew, dim_mess, dim_redund):
	def thre(x):
		if (x > 0):
			return 0
		else:
			return 1

	def mgf2(x):
		return matrix(GF(2),x)

	Lres = []
	for i in xrange(2^dim_redund): # sometimes python seriously sucks 
		Lres.append([])

	for b in range(dim_mess):
		BL = codew.submatrix(row=0,col=b*2^dim_redund,nrows=1,ncols=2^dim_redund)
		BF = BooleanFunction(BL.columns())
		mags = BF.walsh_hadamard_transform()
		for g in xrange(2^dim_redund):
			Lres[g].append(thre(mags[g]))

	Lres = map(mgf2,Lres)

	return Lres

# C = codes.ReedMullerCode(GF(2), 1, #dim_redund)
# get base_redund*m, then use vector encoder, then broadcast, then add bits, should work




# Fast encoder from FMT (or actually, just use the existing Reed-Muller one (tho the fast part is maybe
# less visible). Unfortunately, the Boolean polynomials don't map to the multivariate polynomials over GF(2)...
# oh, but we don't really need to compute an ANF at this point anyways.
# and othw, the unencode method could be used (tho again the performance question could be raised)
# yet we still have a problem of mapping vectors to multivariate polynomials through the canonical embedding
# ah, but that's the "vector" encoding
# would need to add multipoint evaluation for Boolean polynomials using a fast Möbius transform
# also, encoding/decoding Boolean polynomials to and from binary lists. Why is this not even done??
