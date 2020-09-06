try:
    import svgwrite
    import sys
    from Bio import SeqIO
    from os.path import basename
    from numpy import size, array
    import numpy
    import os
except ImportError:
    # if svgwrite is not 'installed' append parent dir of __file__ to sys.path
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.split(os.path.abspath(__file__))[0]+'/..'))
  
import svgwrite
from svgwrite import cm, mm
asblysize=[]

collect=[]
collect2=[]

def calldata():
    header="Assembly","Total Number of Contigs:","Total Length of Contigs:","Total Length of Contigs(excluding N's):","Largest Contig:","Smallest-Contig:","N50:","Contig>=5kb:","Contigs>=5kb length:","Contig>=1kb:","Contigs>=1kb length:","Contig>=500bp:","Contigs>=500bp length:","Contig>=250b:","Contigs>=250b length:","Contig>=100bp:","Contigs>=100bp length:";
    #stats[0][0].add(header)
    collect.append(header)
    #collect2.append(header)
    
   
    for i in range(len(sys.argv)):
                if(sys.argv[i].endswith('.py')):
                            print(sys.argv[i], "is parsing fasta('s)")
                else:
                            filename=sys.argv[i]
                            sizes = [len(rec) for rec in SeqIO.parse(filename, "fasta")]
                            Ncount=0;
                            stats=[]
                            stats2=[]
                            for seq_record in SeqIO.parse(filename, "fasta"):
                                sequence = str(seq_record.seq).upper()
                                Ncount += sequence.count("N")
                                #print("Ncount:",Ncount);
                            all_len=sorted(sizes)
                            asblysize.append(sum(sizes))
                            n2=int(sum(sizes)/2)
                            csum=numpy.cumsum(all_len)
                             
                            # get index for cumsum >= N2
                            csumn2=min(csum[csum >= n2])
                            ind=numpy.where(csum==csumn2)
                            n50 = all_len[int(ind[0])]
                            #print("N50: %s" % n50)
                            file=os.path.splitext(filename)[0]
                            fivekb = [x for x in sizes if x >= 5000]
                            onekb = [x for x in sizes if x >= 1000]
                            onekb2 = [x for x in sizes if (x >= 2500 and x < 5000)] # 1-5kb Newvalue-> (<5000 >2500)
                            halfkb = [x for x in sizes if x >= 500]
                            halfkb2 = [x for x in sizes if (x >= 1000 and x < 2500)] # 500-1kb Newvalue-> (>1000 <2500)
                            quaterkb = [x for x in sizes if x >= 250]
                            quaterkb2 = [x for x in sizes if (x >= 650 and x < 1000)] # 250-500 Newvalue-> (>650 <1000)
                            centbp = [x for x in sizes if x >= 100]
                            centbp2 = [x for x in sizes if (x >= 1 and x < 650)] # <250 Newvalue-> (<650 >1)
                            stats.append(file);                   #Assembly name
                            stats2.append(file)
                            stats.append(place_value(len(sizes)))              #Assembly number of contigs
                            stats2.append(len(sizes))
                            temp=(str(place_value(sum(sizes))) + "bp")
                            stats.append(temp)      #Assembly length bp
                            stats2.append(sum(sizes))
                            temp=str(place_value(sum(sizes)-Ncount)) + "bp (" +  str(round(((Ncount/sum(sizes))*100),3)) + "% N's)"
                            stats.append(temp)    #Assembly length excluding N's
                            stats2.append(temp)
                            stats.append(place_value(max(sizes)))               #Largest contig
                            stats2.append(max(sizes))
                            stats.append(place_value(min(sizes)))                #Smallest contig
                            stats2.append(min(sizes))
                            stats.append(place_value(n50))                     #N50
                            stats2.append(n50)
                            stats.append(place_value(len(fivekb)))                         # Contigs 5kb or >
                            stats2.append(len(fivekb))
                            temp=str(place_value(sum(fivekb))) + "bp (" + str(round(((sum(fivekb)/sum(sizes))*100),2)) + "%)"
                            stats.append(temp)                        # Contigs length 5kb or >
                            stats2.append(sum(fivekb))
                            stats.append(place_value(len(onekb)))                         # Contigs 1kb or >
                            stats2.append(len(onekb2))
                            temp=str(place_value(sum(onekb))) + "bp (" + str(round(((sum(onekb)/sum(sizes))*100),2)) + "%)"
                            stats.append(temp)                        # Contigs length 5kb or >
                            stats2.append(sum(onekb2))
                            stats.append(place_value(len(halfkb)))                         # Contigs 1kb or >
                            stats2.append(len(halfkb2))
                            temp=str(place_value(sum(halfkb))) + "bp (" + str(round(((sum(halfkb)/sum(sizes))*100),2))   + "%)"
                            stats.append(temp)                     # Contigs length 5kb or >
                            stats2.append(sum(halfkb2))
                            stats.append(place_value(len(quaterkb)))                         # Contigs 1kb or >
                            stats2.append(len(quaterkb2))
                            temp=str(place_value((sum(quaterkb)))) + "bp (" + str(round(((sum(quaterkb)/sum(sizes))*100),2))  + "%)"
                            stats.append(temp)                      # Contigs length 5kb or >
                            stats2.append(sum(quaterkb2))
                            stats.append(place_value(len(centbp)))                         # Contigs 1kb or >
                            stats2.append(len(centbp2))
                            temp=str(place_value(sum(centbp))) + "bp (" + str(round(((sum(centbp)/sum(sizes))*100),2))  + "%)"  
                            stats.append(temp)                    # Contigs length 5kb or >
                            stats2.append(sum(centbp2))
                            collect2.append(stats2)
                            collect.append(stats)

def place_value(number): 
    return ("{:,}".format(number))

def transpose(list_in_mat):
    list_out_mat = []
    array_in_mat = array(list_in_mat)
    array_out_mat = array_in_mat.T
    nb_lines = size(array_out_mat, 0)
    for i_line_out in range(0, nb_lines):
        array_out_line = array_out_mat[i_line_out]
        list_out_line = list(array_out_line)
        list_out_mat.append(list_out_line)
    return list_out_mat

def printdata():
    nl=transpose(collect)
    file = open('log.txt','w')
    mx = len(max((sub[0] for sub in nl),key=len))
    
    for row in nl:
        print(" ".join(["{:<{mx}}".format(ele,mx=mx) for ele in row]))
        #file.write(" ".join(["{:<{mx}}".format(ele,mx=mx) for ele in row]))
        for elem in row:
            #print(elem,"UUUUU")
            file.write(elem+"\t")
        file.write("\n")
    file.close()

def createplot(name):
    svg_width=50*cm
    svg_hight=(10*len(sys.argv))*cm
    x=(max(asblysize)/300)
    mid=int(10*len(sys.argv))/len(sys.argv)
    dwg = svgwrite.Drawing(name,(svg_width, svg_hight), debug=True)
    dwg.add(dwg.rect(insert=(0, 0), size=('100%', '100%'), fill='white'))
    lines = dwg.add(dwg.g(stroke_width=5, stroke='green', fill='none'))
    text = dwg.add(dwg.g(font_family="sans-serif", font_size=20, fill='black'))
    text.add(dwg.text("Assembly metrics", insert=(10*cm, (mid-4)*cm), font_size=60, font_weight='bold'))
    shapes = dwg.add(dwg.g(id='shapes', fill='white'))
    curve = dwg.polyline( stroke='green', fill='none', stroke_width=1)
    
    srow=mid-1
    color=["lime","brown","blue","gray","gold","green","maroon","magenta","red","olive"]
    contiglg=["Contigs >5kb","","Contigs 2.5-5kb","","Contigs 1-2.5kb","","Contigs 650bp-1kb","","Contigs <650bp"]
    scol=10
    srow
    #for i in range(1,len(sys.argv)):
    for i in collect2:
        #print(i[2])
        gsize=(i[2]/x)/10
        #for el in i:
         #   print(el,mid,srow)
        scol=10
        #srow=mid
        ecol=(scol+gsize) #37
        erow=srow
        shapes.add(dwg.rect(insert=(scol*cm, srow*cm), size=((ecol-scol)*cm, 3*cm), stroke='black',stroke_width=4 ))
        #lines.add(dwg.line(start=(scol*cm, srow*cm), end=(ecol*cm,erow*cm), stroke='black' ))
        #lines.add(dwg.line(start=(scol*cm, (srow+3)*cm), end=(ecol*cm,(erow+3)*cm), stroke='black' ))
        #lines.add(dwg.line(start=(scol*cm, srow*cm), end=(scol*cm,(erow+3)*cm), stroke='black' ))
        #lines.add(dwg.line(start=(ecol*cm, srow*cm), end=(ecol*cm,(erow+3)*cm), stroke='black' ))
        name=i[0] #Name
        text.add(dwg.text(str(name), insert=((scol-0.5)*cm, (srow)*cm), font_size=24, text_anchor='end',style="font-weight:bold"))
        Length="Length: "+str(place_value(i[2]))
        text.add(dwg.text(str(Length), insert=((scol-0.5)*cm, (srow+0.8)*cm), font_size=24, text_anchor='end'))
        Contigs="Contigs: "+str(place_value(i[1]))
        text.add(dwg.text(str(Contigs), insert=((scol-0.5)*cm, (srow+1.6)*cm), font_size=24, text_anchor='end'))
        n="N50: "+str(place_value(i[6]))
        text.add(dwg.text(str(n), insert=((scol-0.5)*cm, (srow+2.4)*cm), font_size=24, text_anchor='end'))
        Largest="Largest Contig: "+str(place_value(i[4]))
        text.add(dwg.text(str(Largest), insert=((scol-0.5)*cm, (srow+3.2)*cm), font_size=24, text_anchor='end'))
        Smallest="Smallest Contig: "+str(place_value(i[5]))
        text.add(dwg.text(str(Smallest), insert=((scol-0.5)*cm, (srow+4.0)*cm), font_size=24, text_anchor='end'))
        for j in range(8,18,2):
            shpsize=(i[j]/x)
            wid=(shpsize)/10
            shapes.add(dwg.rect(insert=(scol*cm, srow*cm), size=(wid*cm, 30*mm), fill=color[j-8], ))
            #print(x,mid,i[j],scol,ecol,shpsize,wid)
            #srow=srow+0.3
            scol=scol+(wid)
        srow+=mid-3        
    
    #Scale
    #print(srow,svg_hight)
    lines.add(dwg.line(start=(10*cm, (srow-3)*cm), end=(40*cm, (srow-3)*cm), stroke_linecap='round', stroke='black', stroke_width=4))
    for i in range(0,350,50):
        scol=(i/10)+10
        #label=str(i)+" Mbp"
        label=str(round((i*x)/1000000))+"Mbp"
        lines.add(dwg.line(start=(scol*cm, (srow-3)*cm), end=(scol*cm, (srow-2)*cm), stroke_linecap='round', stroke='black', stroke_width=4))
        text.add(dwg.text(str(label), insert=(scol*cm, (srow-1.5)*cm), font_size=20))
    scol=5
    srow=(mid*len(sys.argv))-(len(sys.argv)+2)
    for j in range(8,18,2):
        for i in range(1,10):
            #print(scol,srow,(srow*10)+i)
            lines.add(dwg.line(start=(scol*cm, ((srow*10)+i)*mm), end=((scol+1)*cm, ((srow*10)+i)*mm),stroke=color[j-8]))
        label2=contiglg[j-8]
        text.add(dwg.text(str(label2), insert=((scol+1.1)*cm, ((srow*10)+i)*mm), font_size=24))
        scol=scol+8
        
    dwg.save()

#============
calldata()
printdata()
createplot('AssemblyStats.svg')