### Global tracrRNA predictor
### Prediction method: check the dot-bracket RNA structure pattern, if the gRNA region have multiple same directed bracket, we suppose this is a anti-repeat region
### Use sliding window method to grab tracrRNA sequence, extract all possible tracrRNA.
### tracrRNAPredictionD: for dis-continuous anti-repeat type tracrRNA prediction

import sys,os
import argparse
from Bio import SeqIO
from Bio import pairwise2
import forgi
import forgi.visual.mplotlib as fvm
import matplotlib.pyplot as plt

print("=========================== NOTE ============================")
print("Here, we use UUU as the linker between crRNA & tracrRNA in our secondary structure plot, but in RNAfold file, it will still be LLL.")
print("WARING: Not recommand the condition when len(tracrRNA)+len(gRNA) > 155nt, this will cause forgi error.")

def PostProcessingRNAFoldFile(infile,outfile):
    with open(infile,'r') as fin, open(outfile,'w') as fout:
        finline = fin.readlines()
        fout.write(finline[0].strip()+"\n")
        fout.write(finline[1].replace("LLL","UUU").strip()+"\n")
        fout.write(finline[2].strip().split()[0])
        freeEnergy = finline[2].strip().split()[1]
        DotBracket = finline[2].strip().split()[0]
    return [freeEnergy,DotBracket]

def AcceptableAntiRepeat(inpattern,incomplex,inlength,inthreshold):
    if "LLL" in incomplex:
        LeftFlankSeq, RightFlankSeq = incomplex.split("LLL")
        if len(LeftFlankSeq) == inlength:
            crRNA = RightFlankSeq
        else:
            crRNA = LeftFlankSeq
        crRNAPosition = incomplex.find(crRNA)
        crRNApattern = inpattern[crRNAPosition:crRNAPosition+len(crRNA)]
        ### Now, check the bracket type
        if "(" in crRNApattern and not ")" in crRNApattern:
            ### The bracket direction is "("
            if crRNApattern.count("(") >= inthreshold:
                if "(((" in crRNApattern:
                    return "Acceptable"
                else:
                    return "Ambiguous"
        elif ")" in crRNApattern and not "(" in crRNApattern:
            ### The bracket direction is ")"
            if crRNApattern.count(")") >= inthreshold:
                if ")))" in crRNApattern:
                    return "Acceptable"
                else:
                    return "Ambiguous"
        elif "(" in crRNApattern and ")" in crRNApattern:
            ### Self stem-loop
            AntiRepeatBracketNumber = abs(crRNApattern.count("(")-crRNApattern.count(")"))
            if AntiRepeatBracketNumber >= inthreshold:
                return "Ambiguous"
        else:
            ### No bracket in pattern, crRNA as a free linear structure
            return False
    else:
        return False

def checkEndAR(inseq,intype):
    ARBase3 = ['AA','GA']
    ARBase5 = ['AA','AG']
    if intype == '3':
        far_linkerloc = inseq[:5]
        if 'A' in far_linkerloc:
            for ARBind in ARBase3:
                if ARBind in far_linkerloc:
                    return "PerfectAR"
            Alocation = [i for i, x in enumerate(far_linkerloc) if x == "A"]
            if len(Alocation) == 1:
                FarASide = far_linkerloc[:Alocation[0]]
                if 'A' in FarASide or 'G' in FarASide:
                    return "SaparatedAR"
            else:
                for locA in Alocation:
                    FarASide = far_linkerloc[:locA]
                    if 'A' in FarASide or 'G' in FarASide:
                        return "SaparatedAR"
        else:
            return "NoAR"
    elif intype == '5':
        far_linkerloc = inseq[-5:]
        if 'A' in far_linkerloc:
            for ARBind in ARBase5:
                if ARBind in far_linkerloc:
                    return "PerfectAR"
            Alocation = [i for i, x in enumerate(far_linkerloc) if x == "A"]
            if len(Alocation) == 1:
                FarASide = far_linkerloc[Alocation[0]+1:]
                if 'A' in FarASide or 'G' in FarASide:
                    return "SaparatedAR"
            else:
                for locA in Alocation:
                    FarASide = far_linkerloc[locA+1:]
                    if 'A' in FarASide or 'G' in FarASide:
                        return "SaparatedAR"
    else:
        print("Error: tracrRNA type should be 3 or 5")
        sys.exit(1)

ap = argparse.ArgumentParser()

ap.add_argument("-g","--guide",help="guideRNA sequence", required=True)
ap.add_argument("-t","--tracr",help="putative tracrRNA region, single or multiple sequence fasta file", required=True)
ap.add_argument("-l","--length",help="tracrRNA length", default=130)
ap.add_argument("-o","--outdir",help="output dir", required=True)
ap.add_argument("-k","--keep_tmp",help="keep temp files",action = "store_true")
ap.add_argument("-T","--type",help="tracrRNA type", choices=['3','5'], default='3')
ap.add_argument("-w","--windowsize",help="window size for sliding window", default=50)
ap.add_argument("-b","--bases",help="anti repeat matching base number threshold, default is 7",default=7)

opts = ap.parse_args()

guideRNA = opts.guide
tracrRNAFile = opts.tracr
tracrRNA_length = int(opts.length)
outdir = opts.outdir
keep_tmp = opts.keep_tmp
tracrRNA_type = opts.type
windowSize = opts.windowsize
BaseThreshold = opts.bases

if os.path.isabs(outdir):
    outdir = outdir
else:
    outdir = os.path.abspath(outdir)

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(os.path.join(outdir,"RNAfold/")):
    os.mkdir(os.path.join(outdir,"RNAfold/"))
if not os.path.exists(os.path.join(outdir,"plot/")):
    os.mkdir(os.path.join(outdir,"plot/"))
RNAFoldDir = os.path.join(outdir,"RNAfold/")
PlottingDir = os.path.join(outdir,"plot/")

### Check gRNA string have illegal characters
StandardSequenceBase = ['A','T','G','C','U','a','t','g','c','u']
for bases in guideRNA:
    if not bases in StandardSequenceBase:
        print("Illegal characters in gRNA sequence!")
        sys.exit(1)

fa = open(os.path.join(outdir,"tracrRNAoutput.tab"),'w')
fa.write("### Based on ViennaRNA secondary structure prediction.\n")
fa.write("### hloop: Hairpin loops\n")
fa.write("### mloop: Multi-loops\n")
fa.write("### iloop: Interior loops\n")
fa.write("### floop: Fiveprime regions\n")
fa.write("### tloop: Threeprimes regions\n")
fa.write("###\n")
fa.write("ID\tgRNASeq\ttracrRNASeq\tComplexSeq\tFreeEnergy\thLoop\tmLoop\tStem\tiLoop\tfLoop\ttLoop\tAR\n")

### Now read the tracrRNA region sequence, and use window slide method to find the anti-repeat region
ErrorBox = []
for tracrRecords in SeqIO.parse(tracrRNAFile,'fasta'):
    tracrRegionID = str(tracrRecords.id)
    print("Processing..."+tracrRegionID)
    tracrRegionSeq = str(tracrRecords.seq)
    ### Now, start sliding window, footstep as 1, extract all window sequence
    if len(tracrRegionSeq) < int(tracrRNA_length):
        print("Short tracrRNA, skip..."+tracrRegionID)
    else:
        SlidingWindowLimit = len(tracrRegionSeq) - int(tracrRNA_length) + 1
        for slidingRange in range(0,SlidingWindowLimit):
            tracrRNAID = tracrRegionID+"_Sliding_"+str(slidingRange)
            tracrRNAsequence = tracrRegionSeq[slidingRange:slidingRange+int(tracrRNA_length)]
            ARStatus = checkEndAR(tracrRNAsequence,tracrRNA_type)
            if tracrRNA_type == '3':
                complexRNAseq = tracrRNAsequence + "LLL" + guideRNA
            else:
                complexRNAseq = guideRNA + "LLL" + tracrRNAsequence
            with open(os.path.join(outdir,"RNAfold/",tracrRNAID+".fasta"),'w') as ff:
                ff.write(">"+tracrRNAID+"\n"+complexRNAseq+"\n")
            ff.close()
            RNAFoldCommand = "RNAfold --noPS -i "+os.path.join(outdir,"RNAfold/",tracrRNAID+".fasta")+" > "+os.path.join(outdir,"RNAfold/",tracrRNAID+".ss")
            os.system(RNAFoldCommand)
            tracrRNAfreeEnergy, tracrRNADotBracket = PostProcessingRNAFoldFile(os.path.join(outdir,"RNAfold/",tracrRNAID+".ss"),os.path.join(outdir,"RNAfold/",tracrRNAID+".RNAfold"))
            ComplexRNAStructure = forgi.load_rna(os.path.join(outdir,"RNAfold/",tracrRNAID+".RNAfold"),allow_many = False)
            hloop_num = len(list(ComplexRNAStructure.hloop_iterator()))
            ### Check the hloop number 
            if hloop_num >= 4:
                ### Check the anti-repeat dotbracket
                AntiRepeatStatus = AcceptableAntiRepeat(tracrRNADotBracket,complexRNAseq,tracrRNA_length,BaseThreshold)
                if bool(AntiRepeatStatus):
                    if AntiRepeatStatus == "Acceptable":
                        tracrRNAIDFlagged = tracrRNAID+"_Acceptable"
                        mloop_num = len(list(ComplexRNAStructure.mloop_iterator()))
                        stem_num = len(list(ComplexRNAStructure.stem_iterator()))
                        iloop_num = len(list(ComplexRNAStructure.iloop_iterator()))
                        floop_num = len(list(ComplexRNAStructure.floop_iterator()))
                        tloop_num = len(list(ComplexRNAStructure.tloop_iterator()))
                        plt.figure(figsize=(20,20))
                        fvm.plot_rna(ComplexRNAStructure,text_kwargs={"fontweight":"black"},lighten=0.7,backbone_kwargs={"linewidth":3})
                        try:
                            plt.savefig(os.path.join(outdir,"plot/",tracrRNAIDFlagged+".png"))
                        except:
                            print("tracrRNA complex: "+tracrRNAIDFlagged+" caused a forgi plotting error, this will not output a plot.")
                            print("Please use ViennaRNA RNAfold webpage to draw the RNA structure based on the tracrRNAoutput.tab file.")
                            ErrorBox.append(tracrRNAIDFlagged)
                        fa.write(tracrRNAIDFlagged+"\t"+guideRNA+"\t"+tracrRNAsequence+"\t"+complexRNAseq+"\t"+str(tracrRNAfreeEnergy)+"\t"+str(hloop_num)+"\t"+str(mloop_num)+"\t"+str(stem_num)+"\t"+str(iloop_num)+"\t"+str(floop_num)+"\t"+str(tloop_num)+"\n")
                        plt.close()
                    else:
                        tracrRNAIDFlagged = tracrRNAID+"_Ambiguous"
                        mloop_num = len(list(ComplexRNAStructure.mloop_iterator()))
                        stem_num = len(list(ComplexRNAStructure.stem_iterator()))
                        iloop_num = len(list(ComplexRNAStructure.iloop_iterator()))
                        floop_num = len(list(ComplexRNAStructure.floop_iterator()))
                        tloop_num = len(list(ComplexRNAStructure.tloop_iterator()))
                        plt.figure(figsize=(20,20))
                        fvm.plot_rna(ComplexRNAStructure,text_kwargs={"fontweight":"black"},lighten=0.7,backbone_kwargs={"linewidth":3})
                        try:
                            plt.savefig(os.path.join(outdir,"plot/",tracrRNAIDFlagged+".png"))
                        except:
                            print("tracrRNA complex: "+tracrRNAIDFlagged+" caused a forgi plotting error, this will not output a plot.")
                            print("Please use ViennaRNA RNAfold webpage to draw the RNA structure based on the tracrRNAoutput.tab file.")
                            ErrorBox.append(tracrRNAIDFlagged)
                        fa.write(tracrRNAIDFlagged+"\t"+guideRNA+"\t"+tracrRNAsequence+"\t"+complexRNAseq+"\t"+str(tracrRNAfreeEnergy)+"\t"+str(hloop_num)+"\t"+str(mloop_num)+"\t"+str(stem_num)+"\t"+str(iloop_num)+"\t"+str(floop_num)+"\t"+str(tloop_num)+"\t"+ARStatus+"\n")
                        plt.close()
                else:
                    continue
if len(ErrorBox) != 0:
    print("=================== ERROR ID ================")
    for errorItems in ErrorBox:
        print(errorItems)
