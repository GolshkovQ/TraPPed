### tracrRNA prediction tools
### Predict all putative tracrRNA and draw plot
### tracrRNAPredictionC: for continuous anti-repeat type tracrRNA prediction

from posixpath import isabs
import sys,os
import argparse
from Bio import SeqIO
from Bio import pairwise2
import forgi
import forgi.visual.mplotlib as fvm
import matplotlib.pyplot as plt

print("=========================== NOTE ============================")
print("Here, we use UUU as the linker between crRNA & tracrRNA in our secondary structure plot, but in RNAfold file, it will still be LLL.")


ap = argparse.ArgumentParser()

ap.add_argument("-g","--guide",help="guideRNA sequence", required=True)
ap.add_argument("-t","--tracr",help="putative tracrRNA region, single or multiple sequence fasta file", required=True)
ap.add_argument("-l","--length",help="tracrRNA length", default=130)
ap.add_argument("-o","--outdir",help="output dir", required=True)
ap.add_argument("-k","--keep_tmp",help="keep temp files",action = "store_true")
ap.add_argument("-T","--type",help="tracrRNA type", choices=['3','5'], default='3')
ap.add_argument("-w","--windowsize",help="window size for sliding window", default=50)

opts = ap.parse_args()

guideRNA = opts.guide
tracrRNAFile = opts.tracr
tracrRNA_length = int(opts.length)
outdir = opts.outdir
keep_tmp = opts.keep_tmp
tracrRNA_type = opts.type
windowSize = opts.windowsize

def ReverseComplete(inseq):
    trantab = str.maketrans("ATCGUatcgu","TAGCAtagca")
    return inseq.translate(trantab)[::-1]

def GetContinousAntiRepeatRegion(inID,inSeq,inGuide,inSize,tracrLen):
    ResultBox = {}
    QueryIdentifier = 1
    TempFootStep = 0
    SequencePointerStart = 0
    SequencePointerEnd = SequencePointerStart + int(inSize)
    while SequencePointerEnd <= len(inSeq):
        TargetSequence = inSeq[SequencePointerStart:SequencePointerEnd]
        PairWiseAlignmentRecords = pairwise2.align.globalms(inGuide,TargetSequence,4,0,-10,-0.5)[0]
        AlignmentMarkPattern = str(pairwise2.format_alignment(*PairWiseAlignmentRecords)).split("\n")[1]
        AlignmentScore = str(PairWiseAlignmentRecords[2])
        TempFootStep = GetFootStep(AlignmentMarkPattern,inSize,len(inGuide))
        if TempFootStep[1] != inSize:
            ### This is a putative anti-repeat region
            MatchingAntiRepeatRegionStart = TempFootStep[0]
            MatchingAntiRepeatRegionEnd = TempFootStep[1]
            AntiRepeatRegionStart = SequencePointerStart + MatchingAntiRepeatRegionStart
            AntiRepeatRegionEnd = SequencePointerStart + MatchingAntiRepeatRegionEnd
            if AntiRepeatRegionStart + tracrLen <= len(inSeq):
                AntiRepeatRegionID = inID+"_ARid_"+str(QueryIdentifier)
                QueryIdentifier += 1
                MatchingBases = AlignmentMarkPattern.count("|")
                ResultBox[AntiRepeatRegionID] = [AntiRepeatRegionStart,AlignmentScore,MatchingBases]
            SequencePointerStart += AntiRepeatRegionEnd
            SequencePointerEnd = SequencePointerStart + int(inSize)
        else:
            print("Not found...")
            SequencePointerStart += int(inSize)
            SequencePointerEnd = SequencePointerStart + int(inSize)
    return ResultBox

def GetFootStep(Inpattern,insize,modsize):
    ### If the alignment is acceptable for anti-repeat, let the footstep be the last match length of the alignment
    ### Get first matched position and the last matched position
    if "|||" in Inpattern:
        FirstMatchingPosition = Inpattern.find("|||")
        LastMatchingPosition = -1*Inpattern[::-1].find("|||")
        MatchingLength = len(Inpattern[FirstMatchingPosition:LastMatchingPosition])
        if MatchingLength >= 0.3*modsize and MatchingLength <= 2*modsize:
            return [int(FirstMatchingPosition),int(insize)+int(LastMatchingPosition)]
        else:
            return [0,int(insize)]
    else:
        return [0,int(insize)]

def PostProcessingRNAFoldFile(infile,outfile):
    with open(infile,'r') as fin, open(outfile,'w') as fout:
        finline = fin.readlines()
        fout.write(finline[0].strip()+"\n")
        fout.write(finline[1].replace("LLL","UUU").strip()+"\n")
        fout.write(finline[2].strip().split()[0])
        freeEnergy = finline[2].strip().split()[1]
        '''
            if finline[0] == "." or finline[0] == "(" or finline[0] == ")":
                freeEnergy = finline.strip().split()[1]
                fout.write(finline.split()[0]+"\n")
            else:
                if "LLL" in faline 
                fout.write(finline)
    fin.close()
    fout.close()
    '''
    return freeEnergy

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

ErrorBox = []
fa = open(os.path.join(outdir,"tracrRNAoutput.tab"),'w')
fa.write("### Based on ViennaRNA secondary structure prediction.\n")
fa.write("### hloop: Hairpin loops\n")
fa.write("### mloop: Multi-loops\n")
fa.write("### iloop: Interior loops\n")
fa.write("### floop: Fiveprime regions\n")
fa.write("### tloop: Threeprimes regions\n")
fa.write("###\n")
fa.write("ID\tgRNASeq\ttracrRNASeq\tComplexSeq\tFreeEnergy\tAlignmentScore\tMatches\thLoop\tmLoop\tStem\tiLoop\tfLoop\ttLoop\tAR\n")
### Now read the tracrRNA region sequence, and use window slide method to find the anti-repeat region
for tracrRecords in SeqIO.parse(tracrRNAFile,'fasta'):
    tracrRegionID = str(tracrRecords.id)
    print("Processing..."+tracrRegionID)
    tracrRegionSeq = str(tracrRecords.seq)
    ### Predict continous anti repeat
    AntiRepeatBox = GetContinousAntiRepeatRegion(tracrRegionID,tracrRegionSeq,ReverseComplete(guideRNA),windowSize,tracrRNA_length)
    if len(AntiRepeatBox) == 0:
        print("No anti-repeat region found in "+tracrRegionID)
    else:
        for putAntiID,putAntiValue in AntiRepeatBox.items():
            putAntiPos, putAntiScore, putAntiMatchBase = putAntiValue
            puttracrRNAseq = tracrRegionSeq[putAntiPos:putAntiPos+tracrRNA_length]
            ARStatus = checkEndAR(puttracrRNAseq,tracrRNA_type)
            if tracrRNA_type == "3":
                putativeComplex = puttracrRNAseq+"LLL"+guideRNA
            else:
                putativeComplex = guideRNA+"LLL"+puttracrRNAseq
            with open(os.path.join(outdir,"RNAfold/",putAntiID+".fasta"),'w') as ff:
                ff.write(">"+putAntiID+"\n"+putativeComplex)
            ff.close()
            RNAFoldCommand = "RNAfold --noPS -i "+os.path.join(outdir,"RNAfold/",putAntiID+".fasta")+" > "+os.path.join(outdir,"RNAfold/",putAntiID+".ss")
            os.system(RNAFoldCommand)
            tracrRNAfreeEnergy = PostProcessingRNAFoldFile(os.path.join(outdir,"RNAfold/",putAntiID+".ss"),os.path.join(outdir,"RNAfold/",putAntiID+".RNAfold"))
            ComplexRNAStructure = forgi.load_rna(os.path.join(outdir,"RNAfold/",putAntiID+".RNAfold"),allow_many = False)
            fvm.plot_rna(ComplexRNAStructure,text_kwargs={"fontweight":"black"},lighten=0.7,backbone_kwargs={"linewidth":3})
            hloop_num = len(list(ComplexRNAStructure.hloop_iterator()))
            mloop_num = len(list(ComplexRNAStructure.mloop_iterator()))
            stem_num = len(list(ComplexRNAStructure.stem_iterator()))
            iloop_num = len(list(ComplexRNAStructure.iloop_iterator()))
            floop_num = len(list(ComplexRNAStructure.floop_iterator()))
            tloop_num = len(list(ComplexRNAStructure.tloop_iterator()))
            if stem_num >= 4:
                try:
                    plt.savefig(os.path.join(outdir,"plot/",putAntiID+".png"))
                except:
                    print("tracrRNA complex: "+putAntiID+" caused a forgi plotting error, this will not output a plot.")
                    print("Please use ViennaRNA RNAfold webpage to draw the RNA structure based on the tracrRNAoutput.tab file.")
                    ErrorBox.append(putAntiID)
                fa.write(putAntiID+"\t"+guideRNA+"\t"+puttracrRNAseq+"\t"+putativeComplex+"\t"+str(tracrRNAfreeEnergy)+"\t"+str(putAntiScore)+"\t"+str(putAntiMatchBase)+"\t"+str(hloop_num)+"\t"+str(mloop_num)+"\t"+str(stem_num)+"\t"+str(iloop_num)+"\t"+str(floop_num)+"\t"+str(tloop_num)+"\t"+ARStatus+"\n")
                plt.close()

if len(ErrorBox) != 0:
    print("=================== ERROR ID ================")
    for errorItems in ErrorBox:
        print(errorItems)
