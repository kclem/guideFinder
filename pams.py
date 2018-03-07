import sys
import re
import argparse

def main(arguments):

	parser = argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--exon_file', help="Exon sequence in a file", type=argparse.FileType('r'))
	parser.add_argument('--output_file', help="File to which results are written", type=argparse.FileType('w'))
	parser.add_argument('--exon', help="Exon sequence", type=str)
	parser.add_argument('--stop_codons',help='Comma-separated sequence of stop codons', type=str, default='TGA,TAA,TAG')
	parser.add_argument('--pams',help='Comma-separated sequence of PAMs', type=str, default='NGG,NGA')
	parser.add_argument('--edit_start',help='start index of editing (relative to PAM start)', type=int, default=-17)
	parser.add_argument('--edit_end',help='end index of editing (relative to PAM start) (this base is not included in the editing window)', type=int, default=-12)
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
	print('existinging: ' + str(sorted(existing_stop_codons)))

	fwMatches = find_matches(exon_seq,args.stop_codons,args.pams,args.edit_start,args.edit_end,args.editing_recognition_seq,args.editing_resulting_seq,is_rc=False,locs_to_ignore=existing_stop_codons)
	rvMatches = find_matches(reverse_complement(exon_seq),args.stop_codons,args.pams,args.edit_start,args.edit_end,args.editing_recognition_seq,args.editing_resulting_seq,is_rc=True,locs_to_ignore=existing_stop_codons)

	if args.output_file is not None:
		args.output_file.write('PAM\tediting window start\tediting window end\tedit start\tedit end\tstop codon start\tstop codon frame\tguide start\tguide end\tguide sequence\torientation\n')
	print('PAM\tediting window start\tediting window end\tedit start\tedit end\tstop codon start\tstop codon frame\tguide start\tguide end\tguide sequence\torientation')
	match_count = 0
	fw_match_count = 0
	rv_match_count = 0
	for match in fwMatches:
		(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq) = match
		fw_match_count += 1
		match_count += 1
		print('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\tRV'%(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq))
		if args.output_file is not None:
			args.output_file.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\tFW\n'%(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq))

	for match in rvMatches:
		(match_pam,match_edit_start,match_edit_end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,match_guide_start,match_guide_end,match_guide_seq) = match
		rv_match_count += 1
		match_count += 1
		match_edit_start = len(exon_seq) - match_edit_start
		match_edit_end = len(exon_seq) - match_edit_end
		match_exon_start = len(exon_seq) - match_exon_start
		match_exon_end = len(exon_seq) - match_exon_end
		match_guide_start = len(exon_seq) - match_guide_start
		match_guide_end = len(exon_seq) - match_guide_end


		print('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\tRV'%(match_pam,match_edit_end,match_edit_start,match_exon_end,match_exon_start,stop_codon_start,stop_codon_frame,match_guide_end,match_guide_start,match_guide_seq))
		if args.output_file is not None:
			args.output_file.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\tRV\n'%(match_pam,match_edit_end,match_edit_start,match_exon_end,match_exon_start,stop_codon_start,stop_codon_frame,match_guide_end,match_guide_start,match_guide_seq))

	print('Finished. Found ' + str(match_count) + ' matches (' + str(fw_match_count) +' forward and ' + str(rv_match_count) + ' reverse)')


def find_matches(exon,stop_codons,pams,edit_start,edit_end,editing_recognition_seq,editing_resulting_seq,is_rc=False,locs_to_ignore=[]):
	#object to store matches
	matches = []

	#turn everything to upper
	exon = exon.upper()
	pams = pams.upper()
	stop_codons = stop_codons.upper()

	editing_recognition_seq = editing_recognition_seq.upper()
	editing_resulting_seq = editing_resulting_seq.upper()
	editingRegex = make_regex(editing_recognition_seq)

	#iterate through each PAM
	for pam in pams.split(","):
		pamStr = make_regex(pam)
		pam_locs = re.finditer(pamStr,exon)
		for match in pam_locs:
			start = max(0,min(match.start()+edit_start,len(exon)))
			end = max(0,min(match.start()+edit_end,len(exon)))
			if start > end:
				start,end = end,start
			editing_seq = exon[start:end]
			match_locs = re.finditer(editingRegex,editing_seq)
			for sub_match in match_locs:
				match_exon_start = start + sub_match.start()
				match_exon_end = start + sub_match.end()
				subbed_string = exon[0:match_exon_start] + editing_resulting_seq + exon[match_exon_end:]
				stop_codon_starts = get_stop_codons(subbed_string,stop_codons,is_rc,locs_to_ignore)
				if len(stop_codon_starts) > 0:
					guide_start = max(0,match.start()-20)
					guide_end = min(match.start()+3,len(exon))
					guide_seq = exon[guide_start:guide_end]
					stop_codon_start = stop_codon_starts[0]
					stop_codon_frame = stop_codon_start % 3

					matches.append((pam,start,end,match_exon_start,match_exon_end,stop_codon_start,stop_codon_frame,guide_start,guide_end,guide_seq))
	return matches



def get_stop_codons(seq,stop_codons,is_rc=False,locs_to_ignore=[]):
	stop_codon_arr = stop_codons.split(",")
	ind = 0

	if (is_rc):
		seq = reverse_complement(seq)

	stop_codon_locs = []
	for stop_codon in stop_codon_arr:
		stop_locs = re.finditer(stop_codon,seq)
		for stop_loc in stop_locs:
			if stop_loc.start() in locs_to_ignore:
				continue
			stop_codon_locs.append(stop_loc.start())
	return stop_codon_locs

def contains_stop_codons(seq,stop_codons,is_rc):
	stop_codon_arr = stop_codons.split(",")
	ind = 0
	if (is_rc):
		seq = reverse_complement(seq)
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

	pam_regex_string = "(" + pam_regex_string + ")"

	return pam_regex_string

nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
def reverse_complement(seq):
        return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def reverse(seq):
    return "".join(c for c in seq.upper()[-1::-1])





if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
