import sys
import re
import argparse

def main(arguments):

	parser = argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#			formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--exon_file', help="Exon sequence in a file", type=argparse.FileType('r'))
	parser.add_argument('--output_file', help="File to which results are written", type=argparse.FileType('w'))
	parser.add_argument('--exon', help="Exon sequence", type=str)
	parser.add_argument('--stop_codons',help='Comma-separated sequence of stop codons', type=str, default='TGA,TAA,TAG')
	parser.add_argument('--pams',help='Comma-separated sequence of PAMs', type=str, default='NGG,NGA')
	parser.add_argument('--edit_start',help='start index of editing (relative to PAM start). This base is edited.', type=int, default=-16)
	parser.add_argument('--edit_end',help='end index of editing (relative to PAM start) This base is edited.', type=int, default=-11)
	parser.add_argument('--editing_recognition_seq',help='Recognition sequence for editing', default='C')
	parser.add_argument('--editing_resulting_seq',help='Resulting sequence for editing', default='T')

	args = parser.parse_args(arguments)

	print(args)

	if args.exon is not None:
		exon_seq = args.exon
	elif args.exon_file is not None:
		exon_seq = args.exon_file.readline().rstrip()
	else:
		parser.print_help(sys.stderr)
		sys.stderr.write('exon_file or exon must be specified\n')
		sys.exit(1)

	existing_stop_codons = get_stop_codons(exon_seq,args.stop_codons)

	fwMatches = find_matches(exon_seq,args.stop_codons,args.pams,args.edit_start,args.edit_end,args.editing_recognition_seq,args.editing_resulting_seq,locs_to_ignore=existing_stop_codons)
	rvMatches = find_matches_rc(exon_seq,args.stop_codons,args.pams,args.edit_start,args.edit_end,args.editing_recognition_seq,args.editing_resulting_seq,locs_to_ignore=existing_stop_codons)

	if args.output_file is not None:
		args.output_file.write('PAM\tediting window start\tediting window end\tedit start\tedit end\tstop codon start\tstop codon frame\tguide start\tguide end\tguide sequence\tproduce sequence (stop codon is lowercase)orientation\n')
	print('PAM\tediting window start\tediting window end\tedit start\tedit end\tstop codon start\tstop codon frame\tguide start\tguide end\tguide sequence\tproduct sequence (stop codon is lowercase)\torientation')
	match_count = 0
	fw_match_count = 0
	rv_match_count = 0
	for match in fwMatches:
		(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq,match_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window) = match
		fw_match_count += 1
		match_count += 1
		print('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\tFW'%(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq,match_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window))
		if args.output_file is not None:
			args.output_file.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\tFW\n'%(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq,match_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window))

	for match in rvMatches:
		(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq,match_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window) = match
		rv_match_count += 1
		match_count += 1

		print('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\tRV'%(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq,match_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window))
		if args.output_file is not None:
			args.output_file.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\tRV\n'%(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq,match_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window))

	print('Finished. Found ' + str(match_count) + ' matches (' + str(fw_match_count) +' forward and ' + str(rv_match_count) + ' reverse)')


def find_matches(exon,stop_codons,pams,edit_start,edit_end,editing_recognition_seq,editing_resulting_seq,locs_to_ignore=[]):
	#object to store matches
	matches = []

	exon_upper = exon.upper()
	#turn parameters to upper
	pams = pams.upper()
	stop_codons = stop_codons.upper()

	editing_recognition_seq = editing_recognition_seq.upper()
	editing_resulting_seq = editing_resulting_seq.upper()
	editing_regex = make_regex(editing_recognition_seq)

	#iterate through each PAM
	for pam in pams.split(","):
		pamStr = make_regex(pam)
		pam_locs = re.finditer(pamStr,exon_upper) #search for PAMs in upper
		for pam_match in pam_locs:
			start = max(1,min(pam_match.start()+edit_start,len(exon)))
			end = max(1,min(pam_match.start()+edit_end,len(exon)))
			if start > end:
				start,end = end,start
			editing_seq = exon[start:end+1]
			match_locs = re.finditer(editing_regex,editing_seq)
			for sub_match in match_locs:
				match_exon_start = start + sub_match.start()
				match_exon_end = match_exon_start + len(editing_resulting_seq)
				subbed_string = exon[0:match_exon_start] + editing_resulting_seq + exon[match_exon_end:]
				stop_codon_starts = get_stop_codons(subbed_string,stop_codons,locs_to_ignore)
				if len(stop_codon_starts) > 0:
					guide_start = max(0,pam_match.start()-20)
					guide_end = min(pam_match.start()+len(pam),len(exon))
					guide_seq = exon[guide_start:guide_end]
					stop_codon_start = stop_codon_starts[0][0]
					seq_to_stop_codon = exon[0:stop_codon_start]
					exon_nuc_count = sum(1 for nuc in exon[0:stop_codon_start] if nuc.isupper())
					stop_codon_frame = exon_nuc_count % 3

					subbed_string_list = list(subbed_string)
					for stop_loc in stop_codon_starts[0]:
						subbed_string_list[stop_loc] = subbed_string_list[stop_loc].lower()
					subbed_guide_seq_with_lower = ''.join(subbed_string_list)
					subbed_guide_seq_with_lower = subbed_guide_seq_with_lower[guide_start:guide_end]

					five_prime_base = exon[match_exon_start+1].upper()
					num_cs_in_editing_window = sum(1 for nuc in editing_seq.upper() if nuc == 'C')

					if len(guide_seq) >= 23:
						matches.append((pam,start+1,end+1,match_exon_start+1,match_exon_end,stop_codon_start+1,stop_codon_frame,guide_start+1,guide_end,guide_seq,subbed_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window))
	return matches

def find_matches_rc(exon,stop_codons,pams,edit_start,edit_end,editing_recognition_seq,editing_resulting_seq,locs_to_ignore=[]):
	#object to store matches
	matches = []

	exon_upper = exon.upper()
	#turn parameters to upper
	pams = pams.upper()
	stop_codons = stop_codons.upper()

	editing_recognition_seq_rc = reverse_complement(editing_recognition_seq.upper())
	editing_resulting_seq_rc = reverse_complement(editing_resulting_seq.upper())
	editing_regex_rc = make_regex(editing_recognition_seq_rc)

	#iterate through each PAM
	for pam in pams.split(","):
		pamStr = make_regex(reverse_complement(pam))
		pam_locs = re.finditer(pamStr,exon_upper) #search for PAMs in upper
		for pam_match in pam_locs:
			#get editing window
			start = min(len(exon)+1,pam_match.start()-edit_start+len(pam)-1)
			end = min(len(exon)+1,pam_match.start()-edit_end+len(pam)-1)
			if start > end:
				start,end = end,start
			editing_seq = exon[start:end+1]
			match_locs = re.finditer(editing_regex_rc,editing_seq)
			for sub_match in match_locs:
				match_exon_start = start + sub_match.start()
				match_exon_end = match_exon_start + len(editing_resulting_seq_rc)
				subbed_string = exon[0:match_exon_start] + editing_resulting_seq_rc + exon[match_exon_end:]
				stop_codon_starts = get_stop_codons(subbed_string,stop_codons,locs_to_ignore)
				if len(stop_codon_starts) > 0:
					guide_start = max(0,pam_match.start())
					guide_end = min(pam_match.start()+len(pam)+20,len(exon))
					guide_seq = reverse_complement(exon[guide_start:guide_end])
					stop_codon_start = stop_codon_starts[0][0]
					seq_to_stop_codon = exon[0:stop_codon_start]
					exon_nuc_count = sum(1 for nuc in exon[0:stop_codon_start] if nuc.isupper())
					stop_codon_frame = exon_nuc_count % 3

					subbed_string_list = list(subbed_string)
					for stop_loc in stop_codon_starts[0]:
						subbed_string_list[stop_loc] = subbed_string_list[stop_loc].lower()
					subbed_guide_seq_with_lower = ''.join(subbed_string_list)
					subbed_guide_seq_with_lower = subbed_guide_seq_with_lower[guide_start:guide_end]

					five_prime_base = reverse_complement(exon[match_exon_end-2]).upper()
					editing_seq_rc = reverse_complement(editing_seq)
					num_cs_in_editing_window = sum(1 for nuc in editing_seq_rc.upper() if nuc == 'C')

					if len(guide_seq) >= 23:
						matches.append((pam,start+1,end+1,match_exon_start+1,match_exon_end,stop_codon_start+1,stop_codon_frame,guide_start+1,guide_end,guide_seq,subbed_guide_seq_with_lower,five_prime_base,num_cs_in_editing_window))
	return matches



def get_stop_codons(seq,stop_codons,locs_to_ignore=[]):
	stop_codon_arr = stop_codons.split(",")
	ind = 0

	all_stop_codon_locs = []
	for stop_codon in stop_codon_arr:
		#need to look for stop codons that cross introns
		stop_codon_regex = "(" + ")[a-z]*(".join(list(stop_codon)) + ")"
		stop_codon_locs = re.finditer(stop_codon_regex,seq)
		for stop_codon_loc in stop_codon_locs:
			if stop_codon_loc.start() in locs_to_ignore:
				continue
			loc_tuple = (stop_codon_loc.start(1),stop_codon_loc.start(2),stop_codon_loc.start(3))
			all_stop_codon_locs.append(loc_tuple)
	return all_stop_codon_locs

def contains_stop_codons(seq,stop_codons):
	stop_codon_arr = stop_codons.split(",")
	ind = 0
	while ind < len(seq):
		coding_seq = seq[ind:ind+3]
		if coding_seq in stop_codon_arr:
			return True
		ind += 3
	return False

def make_regex(pamSeq):
	pam_regex_string = pamSeq

	pam_regex_string = pam_regex_string.replace('N','[ATCG]')
	pam_regex_string = pam_regex_string.replace('R','[AG]')
	pam_regex_string = pam_regex_string.replace('Y','[CT]')
	pam_regex_string = pam_regex_string.replace('S','[GC]')
	pam_regex_string = pam_regex_string.replace('W','[AT]')
	pam_regex_string = pam_regex_string.replace('K','[GT]')
	pam_regex_string = pam_regex_string.replace('M','[AC]')
	pam_regex_string = pam_regex_string.replace('B','[CGT]')
	pam_regex_string = pam_regex_string.replace('D','[AGT]')
	pam_regex_string = pam_regex_string.replace('H','[ACT]')
	pam_regex_string = pam_regex_string.replace('V','[ACG]')

	#pam_regex_string = "(?=(" + pam_regex_string + "))"
	pam_regex_string = "(?=(" + pam_regex_string + "))"

	return pam_regex_string

nt_complement=dict({'A':'T',
					'C':'G',
					'G':'C',
					'T':'A',
					'a':'t',
					'c':'g',
					'g':'c',
					't':'a',
					'n':'n',
					'N':'N',
					'R':'Y',
					'Y':'R',
					'S':'W',
					'W':'S',
					'K':'M',
					'M':'K',
					'B':'V',
					'V':'B',
					'H':'D',
					'D':'H',
					'r':'y',
					'y':'r',
					's':'w',
					'w':'s',
					'k':'m',
					'm':'k',
					'b':'v',
					'v':'b',
					'h':'d',
					'd':'h',
					'_':'_',
					'-':'-'})
def reverse_complement(seq):
	return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def reverse(seq):
	return "".join(c for c in seq.upper()[-1::-1])





if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
