from sage.crypto.boolean_function import BooleanFunction

# full encoder; simple but less efficient than doing it in chunks
def PRM_Encoder(mess, dim_redund):
	dim_mess = mess.length()
	if dim_redund > dim_mess:
			raise ValueError

	main_subc	= codes.BinaryReedMullerCode(1, dim_redund)
	subc_E		= main_subc.encoder("EvaluationVector")

	proj_base = codes.random_linear_code(GF(2), dim_mess, dim_redund).generator_matrix().transpose()
	proj_base = block_matrix([zero_matrix(GF(2), dim_mess, 1), proj_base],ncols=2) # add a zero constant coefficient
	proj_mess = mess * proj_base

	main_block 	= matrix(subc_E(proj_mess))
	ones		= ones_matrix(GF(2), 1, 2^dim_redund)
	codew 		= [main_block]*dim_mess
	for i in range(dim_mess):
		if mess[i] == 1:
			codew[i] = codew[i] + ones
	codew = block_matrix(codew,ncols=dim_mess)

	return codew


def FWHT_ListDecoder(codew, dim_mess, dim_redund):
	def vgf2(x):
		return vector(GF(2),x)

	Lres = []
	for i in range(2^dim_redund): # sometimes python seriously sucks
		Lres.append([])

	for b in range(dim_mess):
		BL = codew.submatrix(row=0,col=b*2^dim_redund,nrows=1,ncols=2^dim_redund)
		BF = BooleanFunction(BL.columns())
		mags = BF.walsh_hadamard_transform()
		for g in range(2^dim_redund):
			if mags[g] > 0:
				Lres[g].append(0)
			else:
				Lres[g].append(1)

	return map(vgf2,Lres)


def test(dim_mess, bias):
	mess 		= random_vector(GF(2), dim_mess)
	dim_redund 	= ceil(log(1/(bias^2),2)) #+ 1 for large-ish cases maybe
	print(log(1/(bias^2),2),dim_redund)

	codew = PRM_Encoder(mess, dim_redund)
	print(mess)
	print("*")
	print(codew.ncols())
	print("*")
	error_rate	= 0.5 - bias
	error_pos	= random_matrix(GF(2), 1, codew.ncols(), density=error_rate)
	print(error_rate,vector(error_pos).hamming_weight())
	print("*")
	noisycodew = codew + error_pos

	answ = FWHT_ListDecoder(noisycodew, dim_mess, dim_redund)
	for cm in answ:
		if cm == mess:
			return True

	return False
