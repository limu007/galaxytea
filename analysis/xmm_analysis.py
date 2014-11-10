'''NAME: 	analysis_script
 AUTHOR:	Matteo Guainazzi (RSSD/VILSPA)
      perl version Filip Munz (AsU CZ/ISDC CH)
 VERSIONS: 	1.0 - Feb 13, 2004 - Original Version (from "Sy2_standard_analysis")
               1.1 - Mar 16, 2004 - Support to Timing mode
                                    Fit with PEXRAV continuum
                                    Light curve generation
		1.2 - Oct 31, 2005 - production version of additional features:
				eventlist content checking 
				transfer of results to web interface
				key parameters to SQL database
		1.3	- May 4, 2006 - next milestone passed
				logging corrected 
				observation lists read also from std. input

 PURPOSE:	Script to automatically scan narrow-band features in EPIC spectra
 PARAMETERS    see below

 NOTES

 If number of command line argument < 1 exit
 
'''






## repase TODO

# vyresit err_msg (globalni)
# dokoncit parsing argumentu (taky globalni hash??)

help_text='''
  Usage: analysis_script [mode]
                   [source name] [entry point]
                   [pn threshdold] [MOS threshold] [Obs.ID]
                   [CCF] [Background Dir.]
  [mode]              = mode 	= PN or MOS or BOTH
  [source name]       = name 	= source name
  [Source redshift]   = redsh	= redshift
  [entry point]       = entry	= 1. Six steps to happiness - starting download
				= 2. ODF reprocessing
				= 3. Removing background using GTI
				= 4. Generating spectra (+ image) + response matrices
				= 5. Background + lightcurves, rebinning spectrum
				= 6. Spectral analysis
  [bkg threshold]     = lim_pn
  			lim_mos	= background rejection threshold
  [Obs.ID]            = obsid	= observation number
  [CCF]               = ccf_dir	= CCF location
  [Background Dir]    = bkg_dir = Directory where Read & Ponman background files are located

  other parameters:
                        proc = run really ?
                        loud = print more output
                        clean = delete temp. files on exit
                        pict = make PNG pictures
                        sql = print SQL commands
						anal = alternative analysis of spectra in QDP files
						check = level of checking for intermediate results
						
						rasc, decl = position of the searched object
'''

msg={} #'err':''} #report system

source=""

redshift=0
(ra_pnt,dec_pnt)=(0,0) # pointing direction
(ra_obj,dec_obj)=(0,0) # catalog object location
pa_pnt=-1
nh=0
goodtime=0
global evtscnt
evtscnt=0
evtsmin=1000
lim_rate={'pn':1.,'mos':1.}

global work_dir
work_dir=""
arcaddr="xsa.vilspa.esa.es"
regime="Imaging"
outip="limu\@zeus.asu.cas.cz";
outaddr=outip+":/var/www/xmm/res"
#$outaddr="sirrah.troja.mff.cuni.cz:www/agn";

jdkhome="/usr/lib/jvm/java";
srchome="/home/limu/bin";
lheahome="/opt/share/heasoft-6.15";
calibhome="/opt/xmm/sas/ccf";

# data repository
archome="/var/data/xmm"; # where downloaded tar files are stored
#archome="/home/munz/data"; # where downloaded tar files are stored
runhome="/home/limu/Lab/xmm/runs"; # where runs are processed
arcresp="/home/limu/Lab/xmm/resp"; # where RMFs and ARFs are backed up
#$obs_table="/home/limu/cat/list_AGN_sorted.txt";
obs_table="/home/limu/Lab/xmm/list_final.txt";
ver_dir="/home/limu/Lab/xmm/veron/";
ccf_directory=calibhome+"/constituents"

modnames={"po":"power-law","pozg":"power-law with 1 gaussian",
    "pozgzg":"power-law with 2 gaussians","pexrav":"power-law with exponential cut-off"}
    
bore_sky_pos={'pn':"25900,25900","mos":"25900,25900"}

ccf_directory=calibhome+"/constituents"
bkg_directory=None

modtxt='both'
global is_inp,oufile
global loghand

is_inp=0
oufile=None
loghand=None

import pyfits
import os
import re
from string import upper as uc
from string import lower as lc
from glob import glob
env=os.environ
env["JDK_HOME"]="/usr/local/j2sdk/"
env["DATA_REP"]=archome
env["SAS_CCFPATH"]=ccf_directory
if (not bkg_directory) :
    bkg_directory = calibhome+"/background"
if ("XMM_PROC" in env) :
        gproc=env["XMM_PROC"]
if ("XMM_PICT" in env) :
        gproc=env["XMM_PICT"]

def get_patt(text,patt,idx=0):
    res=re.search(patt,text)
    if res:
        if idx==None: return res.groups()
        return res.groups()[idx]

root_dir = os.getcwd().strip()
rep=get_patt(root_dir,"\/[^\/\d]+?(\d{9,})")
if (rep) :
    if not obs_number: obs_number=rep
    #root_dir=$`

        
def get_keys(extname,ilist,ext=0):
    status = 0
    extname=extname.replace(":(\w+)","\[\1\]");
    try:
        fptr=pyfits.open(extname)
    except:
        print "Cannot open FITS file ",extname
        fptr=0;status=0
        return []
    olist=[ftpr[ext].header.get(l,'') for l in ilist];
    return olist;

def get_column(extname,colnum,ext=0):
    try:
        fptr=pyfits.open(extname)
    except:
        print "Cannot open FITS file ",extname
        return
    return fptr[ext].data.column(colnum)
    
def get_nh(ra_obj,dec_obj):
    val=0
    rep=os.popen("nh map=%s/refdata/h1_nh.fits equinox=2000 ra=%f dec=%f"%(lheahome,ra_obj,dec_obj))
    out=rep.search("average.+?([\+\-\d\.E]+)\s*$")
    if out:
        val=val.groups()[0]
    return val/1e22

def ifind(text,patt):
    return lc(text).find(lc(patt))>=0
    
def print_head(txt):
    print txt
    
def print_sys(cmd,loud=1):
    #log and execute
    os.system(cmd)
    
def get_sys(cmd):
    return os.popen(cmd).read()
    
####################################################

def init_obs(obs_number):
    global loghand
    if not obs_number:
        print "Problem in "+base_dir+"\n"
        print_help()
        return -1
    
    base_dir=root_dir+"/"+obs_number
    if not os.isdir(base_dir):
        os.mkdir(basedir)
    print "Base directory: $base_dir\n"; 
# 	chomp $base_dir;
    env["SAS_ODF"]=base_dir+"/ODF"
    
    if (not gproc and gloud<2): gloud=2
    if not source:
        OBSL=open(obs_table)
        for li in OBSL.readlines():
            if not li.find("\| "+obs_number+" \|"):
                continue
            data=[a.strip() for a in li.split("|")]
            rpos=[i for i in range(len(data)) if data[i]!=""]
            if len(rpos)>0:
                data=data[rpos[0]:]
            #print "DATA: @data\n";
            source,redshift,coord=tuple(data[:3])
            obstime,obsexpos=data[9],data[11]
            from numpy import array
            amul=array([1,1/60.,1/60.**2])
            m1=re.search(coord,"(\d+):(\d+):(\d+)")
            if (m1):
                ra_obj=(array(map(float,m1.groups()))*amul).sum()*15
                m2=re.search(coord[m1.end():],"(\d+):(\d+):(\d+)")
                if (m2) :
                    dec_obj=(array(map(float,m2.groups()))*amul).sum()
                    if coord.find("-")>=0: dec_obj=-dec_obj

            ra_obj=round(ra_obj,3)
            dec_obj=round(dec_obj,3)
            obstype=" "
            if ifind(data[13],"pn"): obstype+="pn " 
            if ifind(data[13],"mos"): obstype+="mos " 
            
            source=source.rstrip()#replace(" ","")
            print "Found obs ["+obstype+"] at ("+ra_pnt, dec_pnt+") : from "+obstime+" exposure "+obsexpos+"\n"
        OBSL.close()
        if not source:
            print "Cannot find observation no. "+obs_number+".. in "+obs_table+"\n"
            return -2
    corrsrcname='''
    $source=~s/[\.\,]//g;
    $source=~s/\+/_/g;
    $source=~s/ +/_/g;
    $source=~s/_$//g;
    $source=~s/([_-])[_-]+/$1/g;
    #maximum length
    $source=~s/^(.{10,14})_.+/$1/;
    $source=~s/MARK/Mkn/i;
    $source=~s/Mrk/Mkn/i;
    '''
    
    #if ra_obj.find(":")>0:
        #$ra_obj=($1+($2+$3/60.)/60.)*15 if $ra_obj=~/(\d+):(\d+):(\d+)/;
        #$dec_obj=$1+($2+$3/60.)/60. if $dec_obj=~/(\d+):(\d+):(\d+)/;
        #$ra_obj=~s/(\.\d{3})\d+/$1/;
        #$dec_obj=~s/(\.\d{3})\d+/$1/;
    
    if not ra_pnt:
        ra_pnt=ra_obj
        dec_pnt=dec_obj
    if (gloud>0):
        print "Mode "+modlst+" ["+gproc+"] on "+obs_number+" : "+source+" [%.3f,%.3f] \n"%(ra_obj , dec_obj);
        print "\t redshift "+redshift+" bkg "+lim_rate[mod]+" \n";

    if (glog):
        loghand=open(base_dir+"/"+source+"_out.log","w")
        loghand.write("started "+today)


def gen_rmf_arf():
    bkpset=source+"_"+mod+".rmf"
    matset=work_dir+"/"+bkpset
    bkpset=arcresp+"/"+bkpset
    speset="spectrumset="+work_dir+"/"+source+"_"+mod+".pi";
    if os.path.exists(matset) :
        print "Matrix "+matset+" already exists..\n"
    elif os.path.exists(bkpset) :
        print "Linking "+matset+" to "+bkpset+"..\n";
        os.symlink(bkpset,matset)
    else :
        print_head("Generating $mod redistribution matrix")
        print_sys("nice rmfgen rmfset="+matset+" "+speset+"; mv "+matset+" "+bkpset)
        is_ok=1
        if (not os.path.exists(bkpset) and gproc):
            is_ok=0
            msg["err"]="cannot generate "+matset
        else :
            os.symlink(bkpset,matset)

    bkpset=source+"_"+mod+".amf"
    matset=work_dir+"/"+bkpset
    bkpset=arcresp+"/"+bkpset
    if os.path.exists(matset) :
        print "Matrix "+matset+" already exists..\n"
    elif os.path.exists(bkpset) :
        print "Linking "+matset+" to "+bkpset+"..\n";
        os.symlink(bkpset,matset)
    else :
        print_head("Generating $mod arf table")
        print_sys("nice arfgen arfset="+matset+" "+speset+"; mv "+matset+" "+bkpset)
        is_ok=1
        if (not os.path.exists(bkpset) and gproc) :
            is_ok=0
            msg["err"]="cannot generate "+matset
        else:
            os.symlink(bkpset,matset) #move to response directory and link

def pharbn(instr,phafile,bkgfile,resmat,nchan,ncounts,syserr):

    rootpha=get_patt(phafile,"(.+)\.")
    resmat_root=get_patt(resmat,"(.+)\.")
    ancrfile = resmat_root+'.arf'
     #nchan  'Number of channels for resolution element >'
    # ncounts  'Number of minimum counts per bin >'

    if imod.find("PN")>=0:
        res = 2.26 #'Instrumental resolution at Reference energy keV >'
        eref = 6.0 # 'Reference energy (keV) '
        expo = 0.50 # 'Exponent of resolution vs energy law >'
        flchan = "0-4095" #
        eoffset = 0.10222 # 'Energy-indepent resolution component (XMM-case)'
    else:
        res = 6.31
        eref = 1.0
        expo = 0.45702
        flchan = "0-799"
        eoffset = 0.2079893
    
    
    emin=get_column(resmat,2,2)
    emax=emin[1:]
    emin=emin[:-1]
    # simple and dirty
    nchan=len(emin)
    ener=(emin+emax)/2.
    
    (ttype,npha,srcback)=get_keys(phafile,['TTYPE','NAXIS2','BACKSCAL'],1)
    srcbins=get_column(phafile,2,1)
    (ttype,npha,bkgback) =get_keys(bkgfile,['TTYPE','NAXIS2','BACKSCAL'],1)
    bkgbins=get_column(bkgfile,2,1)
    srcbins-=bkgbins*srcback/bkgback
    
    # resampling channels in new exponential deal
    enew=[emin[0]]
    c1=res/100.*eref;
    j=0
    enext=0
    for i in range(nchan):
        if (enext==0 or emax[i]<enext):
            if enext>0:
                enew.append(emax[i])
            enext=enew[j]+c1*exp((eoffset+enew[j]/eref)*expo)/nchan
            j+=1
    nchanew=len(enew);
    
    #newbin=[srcbin[lev==i] for i in range(nchanew)]
    newbin=[srcbin[(ener>enew[i])*(ener<enew[i+1])] for i in range(nchanew-1)]
    #j=0;
    #for (my i=0;i<npha;i+=1) {
        #if (ener[i]>enew[j+1]) {
            #last if j>nchanew-1;
            #j+=1;
            #newbin[j]=0.;
        #}
        #newbin[j]+=srcbin[i];
    #}
    efin=[]
    sumbin=0
    for i in range(nchanew):
        sumbin+=enew[i]
        if sumbin>ncount:
            sumbin=0
            efin.append(enew[i])
    efin.append(enew[-1])
    
    
    #my grouping = rootpha.'.grouping';
    outfile=open(rootpha+'.grouping',"w")# or die "cannot write rebinned bands!";
    j=0
    lincnt=0;
    ipas=ireb1=ireb2=1;
    for i in range(nchan):
        if (ener[i]>efin[j] and ener[i]<=efin[j+1]) :
            if ipas==1:
                ireb=i
            ipas+=1
        elif (ener[i]>efin[j+1] and ener[i]<=efin[j+1]) :
            ireb2=i-1
            nreb=ireb2-ireb1+1
            if nreb>1:
                outfile.write("\t".join([ireb1,ireb2,nreb]))
                lincnt+=1
            j+=1

    ireb2=i
    nreb=ireb2-ireb1+1
    outfile.write("\t".join([ireb1,ireb2,nreb]))
    lincnt+=1
    outfile.close()
    print "Final "+lincnt+" channels\n"
    
    outphafile = rootpha+'_r.pi';
    extra="";
    for di in open(rootpha+'.cmd').readlines():
        extra+=' '+di.strip();
    
    comm="grppha infile=%{phafile} outfile=%{outphafile} comm=\"reset grouping "
    comm+="& group %{grouping} & chkey RESPFILE %{work_dir}/%{resmat}"
    comm+=" & chkey ANCRFILE %{work_dir}/%{ancrfile} & chkey BACKFILE %{work_dir}/%{bkgfile}"
    comm+=" & systematics %{flchan} %{syserr} & %{extra} & exit\" clobber=yes"
    print comm%locals()
    
    test="data 1 "+outphafile+"\ncpd /xt"
    test+="setplot rebin 3 4096\nsetplot channel\n"
    test+="query yes\nplot counts\nsetplot energy"
    try:
        ftest=open(resmat_root+".xcm","w")
        ftest.write(test+"\n")
        ftest.close()
    except:
        print "failed writing "+resmat_root+".xcm"

def bin_spectrum(name,*pars):
    print_head("Spectral grouping")
    #my $pars=join " ",@_;
    emod='E'+uc(mod);
    if mod.find("mos")>=0: 
        emod+='1'
    is_ok=1
    os.chdir(work_dir)
    comp=srchome+"/pharbn "+emod+" \""+name+".pi\" \""+name+"_bkg.pi\" "+name+".rmf "+" ".join(pars)
    print_sys(comp)
    if is_ok and not os.path.exists(work_dir+"/"+name+"_r.pi"):
        is_ok=0
        msg["err"]="rebinning didn't work"

############### XSPEC part ####

def read_inp(source):
    inpfile=work_dir+"/"+source+".inp";
    infile=open(inpfile)# or print "Cannot open input data $inpfile\n";
    for li in infile.readlines():
        imod=get_patt(li,"(\w+) PROPERTIES")
        rep=get_patt(li,"coordinates: RA=([\d\.\+\-]+), Dec=([\d\.\+\-]+), PA=([\d\.\+\-]+)")
        if rep:
            ra_pnt,dec_pnt,pa_pnt=rep
        elif (imod) :
            rep=get_patt(li,"region coordinates: \(.([\d\.\+\-]+,[\d\.\+\-]+).,.([\d\.\+\-]+).\)")
            if (rep) :
                bore_sky_pos[imod],bore_sky_rad=rep
                print "re-using position "+bore_sky_pos[imod]+" ("+bore_sky_rad+")\n";
            elif ifind(li,"FILTER") :
                ifilter[imod]=get_patt(li,"FILTER +(\w+)")
            elif ifind(li,"SUBMODE") :
                submode[imod]=get_patt(li,"SUBMODE +(\w+)")
    infile.close()

def write_inp(source,imod):
    global is_inp,oufile
    oufile=open(work_dir+"/"+source+".inp","w")# or print "Cannot open inpfile\n")
    oufile.write("SOURCE PROPERTIES\n" )
    oufile.write("-----------------\n" )
    oufile.write("Source name: "+source+" \n")
    oufile.write("Median observation coordinates: RA="+ra_pnt+", Dec="+dec_pnt+", PA="+pa_pnt+" \n")
    oufile.write("Galactic N_H (10^22 cm^-2): "+nh+" \n")
    oufile.write("Source redshift: "+redshift+" \n\n")
    
#	for (my imod in ["pn","mos"]) {
    if (imod) :
        oufile.write(imod+"PROPERTIES\n")
        oufile.write("-------------\n")
        oufile.write(imod+"instrument configuration:\nFILTER "+ifilter[imod]+" \nSUBMODE "+submode[imod]+"\n")
        oufile.write(imod+"extraction region coordinates: (\""+bore_sky_pos[imod]+"\",\""+bore_sky_rad+"\")\n")
        oufile.write(imod+"background threshold (cts/s): "+lim_rate[imod]+"\n")
    is_inp=1
    #oufile.close()

# continuum norm, , , , peak width
parval=(1.0,	2.0,	0.5,	2.0,	1.0,	1.0E-2);
pardif=(1.0E-2,1.0E-2,	-1.0E-2,-1.0E-2,1.0E-2,	1.0E-3);
parlow=(0.0,	-3.0,	-100.0,	-100.0,	0.0,	0.0);
parhigh=(1.0E6,10.0,	1.0E10,	1.0E10,	1.0E24,	20.0);

def print_par(fval,fdif,flow=0,fhigh=1E24,fmin=None,fmax=None):
    # formatted output of XSPEC parameters
    if not fmin: fmin=flow
    if not fmax: fmax=fhigh
    
    sformat="%12.4f %12.4e %12.4f %12.4f %12.4e %12.4e\n"
    if fval!=0 and abs(fval<1E-3):
        sformat="%12.4e %12.4e %12.4f %12.4f %12.4e %12.4e\n"
    xcmfile.write(sformat,(fval,fdif,fmin,flow,fhigh,fmax))

spec_int_wide="3.0 12.0"
spec_int_narrow="4.0 7.0"

def get_rates(source):
    nbins=0
    print_head("I am calculating count rates");
    os.chdir(work_dir)
    print_sys("xspec - "+source+"_load.xcm",0)
    filname=source+"_load.log";
    logfile=open(work_dir+"/"+filname)
    log=logfile.read()
    countrate=get_patt(log,"observed count rate +([\d\.]+)")
    #$countrate+=$1 if /(\+\/\-[\d\.eE\+\-]+)/;
    slope==get_patt(log,"PhoIndex +([\d\.eE\+\-]+)")
    norm=get_patt(log,"norm +([\d\.eE\+\-]+)") 
    counts=get_patt(log,"Source file counts +: +(\d+)")
    nbins=get_patt(log,"using +(\d+) +PHA bins")
    logfile.close()
    if (gloud>1) :
        print "normalisation ",norm,"exp ",slope," gives ",counts," counts [",nbins," bins]\n"
    exposure[mod]=int(get_keys(work_dir+"/"+source+"_"+mod+".pi",["ontime"],2))
    if is_inp :
        oufile.write(mod+"exposure time: "+exposure[mod]+" (s) \n")
        oufile.write(mod+" count rate: "+countrate+"\n")

####################################################

def get_data(number): # data retrieval using JAVA
    print_head("Retrieveing the ODF via AIO")
    if os.path.exists(archome+"/"+number+".tar.gz"):
        print "File ",number,".tar.gz already downloaded - extracting"
    else :
        print "File ",number,".tar.gz not present - downloading"
        comp=srchome+"/aioclient "
        if jdkhome:
            comp+="-jdkhome="+jdkhome
        comp+="-S "+arcaddr+" -P 2002 -L \'GET obsno="+number+" level=ODF\' ;"
        comp+="mv "+number+".tar.gz "+archome+"/";
        print_sys(comp,0)
    from glob import glob
    arclst=glob(archome+"/"+number+"*.tar.gz")
    if len(arclst)==0:
        is_ok=0
        msg["err"]="no downloaded file found!"
        return
    if not os.path.isdir(base_dir+"/ODF"):
        os.mkdir(base_dir+"/ODF") 
#	$comp="cd $base_dir/ODF";
    os.chdir(base_dir+"/ODF")
    print_sys("tar xzvf "+arclst[0])
    arclst=glob("*.TAR")
#	@arclst=glob $base_dir."/ODF/*.TAR";
    if len(arclst)==0:
        is_ok=0
        msg["err"]="no archive included in the downloaded file!"
        return

    print_sys("tar xvf "+arclst[0])
    os.unlink(arclst[0])

def proc_data(pmod):
    logfile=source+"_"+comp+"_"+processing_date+".log"
    os.chdir(base_dir+"/"+pmod)
    print_sys("nice "+(pmod.find("pn")>=0 and "ep" or "em")+"proc")

def find_data(mod,step=1,nmin=10): 
    '''
    what came out of the data pipeline
    '''
    global evtscnt
    if ((step>=4 and step<7) and gproc<3):
        filelst=glob(base_dir+"/spectra/*"+mod+".evt")
        if len(filelst)>0:
            regime="Final";
            outfile=filelst[0];
            naxis2,telapse,livetime = get_keys(outfile,["naxis2","telapse","livetime"],ext=1);
            evtscnt=naxis2;
            goodtime=livetime;
            print "Using data from ",outfile," [step ",step," proc ",gproc,"]"
            return outfile
    

    filelst=glob(base_dir+"/"+mod+"/*_ImagingEvts.ds")
    if len(filelst)>0 : #some exist
        regime="Imaging";
    else :
        regime="Timing";
        filelst=glob(base_dir+"/"+mod+"/*_TimingEvts.ds")

    exe,exreg=(env["HEADAS"]+"/bin/ftstat",srchome+"/extract_"+mod+".reg");
    if (not os.path.exists(exe)):
        exe=''
    if regime.find("Imaging")>=0:
        if os.path.exists(exreg): 
            exreg="[regfilter(\'"+exreg+"\',X/20,Y/20) ]"
        else:
            exreg=''
    else:
        exreg=''
    exbin="[bin (X,Y)=::20 ]"
    #$exreg='';
    evtscnt=0;sumtime=0;
    tmpevts={}
    for file in (sorted(filelst)):
        if file.find("E"+uc(mod))==-1: continue
        if file.find("_WORK")>=0: continue
        (naxis2,telapse,livetime,submod) = get_keys(file,["naxis2","telapse","livetime","submode"],1)
        if naxis2<=nmin:
            print file," contains no events!\n"
            continue
        if file.find("_MERG")>=0:
            filelst=(file,)
            evtscnt=naxis2
            break
        evtscnt+=int(naxis2)
        sumtime+=int(livetime)
        tmpevts[file]=naxis2

    filecnt =len(filelst)
    if filecnt>1:
        filelst = tmpevts.keys()
        filecnt = len(filelst)
    
    if filecnt==0:
        #$is_ok=0;
        msg["err"]="No usable "+mod+" files found..\n"
        return
    
    if evtscnt<evtsmin:
        #$is_ok=0;
        msg["err"]="Too few events ("+evtscnt+") found..\n"
        return
    
    print "Found "+filecnt+" "+mod+" event file(s) with "+evtscnt+"events..\n"
    goodtime=sumtime
    is_ok=1
    if filecnt>1:
        workevts=tmpevts[filelst[0]]
        patt=re.compile("(.+)_E"+mod+"(\d*)",re.I)
        res=patt.finditer(infile).next()
        if res:
            d1,d2=res.groups()
            print "Merging ",mod," files .."#;STDOUT->autoflush(1);
            prefix=d1+"_E"+uc(mod)
            if (d2): 
                prefix+='0'
            workfile=prefix+"_WORK_"+regime+"Evts.ds"
            outfile=prefix+"_MERGE_"+regime+"Evts.ds"
        
        i=0
        for actfile in (sorted(filelst)):
            comp="merge set1=\""+infile+"\" set2=\""+actfile+"\" outset=\""+workfile+"\";"
            comp+="mv "+workfile+" "+outfile
            (newevts,submod)=get_keys(actfile,["naxis2","submode"],1);
            for a in re.findall("\/"+mod+"\/(.+)",actfile):
                print "including "+a+" (",submod," - ",newevts,")\t"
            if (exe and exreg and gcheck>0):
                maxpic=0;
                com_check=exe+' "'+actfile+'[EVENTS]'+exreg+exbin+'"'
                res_check=os.popen(com_check).read()
                if ganal>2:
                    print com_check,"\n",res_check;
                else :
                    #while ($res_check=~s/mean:.+?([\d\.]+)//) {$pos+=$1.", ";}
                    res=re.search("maximum.+?(\d+)",res_check)
                    if res: 
                        maxpic=res.groups()[0]
                    res=re.search("sum.+?(\d+)",res_check)
                    if (res!=None and newevts):
                        resevts=round(100*(1-res.groups()[0]/newevts),2)
                    else: 
                        resevts=0
                    print resevts,"% passed [",maxpic,"]";
                
            
            infile=outfile
            print_sys(comp)
            if not os.path.exists(outfile):
                is_ok=0
                msg["err"]="Merging of "+mod+" files failed!"
                break
            elif gcheck>1:
                i+=1
                print i," - ",newevts
                if (newevts!=workevts+tmpevts[actfile]):
                    msg["err"]="Number of events changed during merging : "+newevts+" != "+workevts+" + "+tmpevts[actfile]
                    if os.path.exists(outfile):
                        (newevts,)=get_keys(outfile,["naxis2"],1) #checking again
                    newevts-=workevts;
                    (workevts)=get_keys(actfile,["naxis2"],1) 
                    newevts-=workevts;
                    if (abs(newevts)>100):
                        is_ok=0;
                        break
                else:
                    print ".OK.";#STDOUT->autoflush(1);
                    workevts=newevts
                
            
        
        print 
    else:
        outfile=filelst[0]

    if is_ok and regime.find("Imaging")>=0 and ra_pnt==ra_obj:
        var_list=get_keys(outfile,["ra_pnt","dec_pnt","pa_pnt"],1)
        (ra_pnt,dec_pnt,pa_pnt)=[round(float(var),3) for var in (var_list)]
        print "Pointing position %5.3f , %5.3f [target offset %5.2f , %5.2f arcmin]\n"%(ra_pnt,dec_pnt,(ra_pnt-ra_obj)*60.,(dec_pnt-dec_obj)*60.)
        # problem with use strict refs..
        #foreach my $var (@var_list) {$$var=$$var*1.00;$$var=~s/(\d+\.\d{3})\d+/$1/ ;}
    
    return outfile
    
def clean_data(*argv):
    rep=os.popen("du -s "+base_dir).read()
    res=re.search("^ *(\d+)",rep)
    if res:
        size=res.group()[0]
    comp="cd "+base_dir+";rm -rf ODF ";
    dellst=[]
    for lmod in (argv):
        dellst.extend(glob(base_dir+"/E1*"+uc(substr(lmod,0,1))+"?.fits"))
        dellst.extend(glob(base_dir+"/light_curves/*"+lc(lmod)+"_bs_*_dt10.lcu"))
        comp+=lmod+" "
        if (gsql>0) :
            file=find_data(lmod,step)
            quest="update xmm_anal set evtraw="+evtscnt+" where obsid=\""+obs_number+"\" and mode=\""+mod+"\";\n";
            # replaced by the following (assuming MySQL > 5.0)
            quest="insert into xmm_anal (obsid,mode,evtraw) values ("+','.join(obs_number,mod,evtscnt)+")"
            quest+=" ON DUPLICATE KEY UPDATE set evtraw="+evtscnt
            sqlhand.write(quest)
    dellst=glob(base_dir+"/*SUM.SAS")
    
    print_sys(comp)
    for d in delst:
        os.unlink(d) 
    rep=os.popen("du -s "+base_dir).read()
    res=re.search("^ *(\d+)",rep)
    if res:
        print "Reduced from "+size+" kbytes to "+res.group()[0]
        size=res.group()[0]
    
#############################################################

def gen_cif(name,redo=0):
    env["SAS_CCF"]=base_dir+"/"+name
    if not os.path.exists(env["SAS_CCF"]):
        env["SAS_CCF"]=base_dir+"/ODF/"+name 
    is_cif=0
    if (redo<0 and not os.path.exists(env["SAS_CCF"])):
        return
    if (redo>0 or not os.path.exist(env["SAS_CCF"])):
        print_head("Generating CIF")
        comp="cifbuild withccfpath=yes fullpath=yes "
        comp+="ccfpath="+env["SAS_CCFPATH"]+" withobservationdate=yes"
        print_sys(comp)
        env["SAS_CCF"]=base_dir+"/"+name
        env["SAS_CCF"]=base_dir+"/ODF/"+name 
        if not os.path.exists(env["SAS_CCF"]):
            env["SAS_CCF"]=base_dir+"/ODF/"+name 
    is_cif=1


def ingest_odf(proc=0):
    sumlst=glob(base_dir+"/*SUM.SAS*")
    if len(sumlst)==0:
        sumlst=glob(base_dir+"/ODF/*SUM.SAS*")
    if (proc or len(sumlst)==0):
        if not os.path.exists(base_dir+"/ODF/"):
            return
        env["SAS_ODF"]=base_dir+"/ODF/"
        if gproc and len(sumlst)>0:
           os.unlink(sumlst[0])
        os.chdir(base_dir)
        print_sys("odfingest")
        sumlst=glob(base_dir+"/*SUM.SAS*")
    
    #	print "SUMLST:@sumlst\n";
    env["SAS_ODF"]=sumlst[0]
    if gloud>2: 
        print "odf summary set to ",env["SAS_ODF"],"\n"
    res=re.search("[^\/]+\.SAS",env["SAS_ODF"])
    if res:
        return res.group()[0]

################ evselect routines ###########################

binsize=20
timebinsize=10
specbinsize=5
imgfrac=0.001

def gen_pileup(tabname,extexpr):

    gtifile=source+"_"+mod+"_bkg.gti"
    print_head("Pile-up estimation")
    #my $comp="cd $work_dir/;";
    comp=""
    expr=""

    if (not os.path.exists(work_dir+"/"+source+"_"+mod+".evt")) or gproc>1:
        comp="evselect table="+tabname+" "
        comp+="withfilteredset=yes filteredset=\""+source+"_"+mod+".evt\" keepfilteroutput=yes "
        if os.path.exists(work_dir+"/"+gtifile):
            expr="gti("+gtifile+",TIME)"
        if extexpr: 
            if expr:
                expr+="&&"+extexpr
            else:
                expr=extexpr
        comp+="expression=\""+expr+"\""
        print_sys(comp)

    (naxis2,ra_nom,dec_nom) = get_keys(work_dir+"/"+source+"_"+mod+".evt",["naxis2","ra_nom","dec_nom"],1)
    if evtscnt>0:
        print "Extracted %i evts = %3.1f%% [%s]\n"%(naxis2,naxis2/evtscnt*100.,extexpr)
    else:
        if naxis2>0:
            print "No events!\n"
        else:
            print "Events of no source!!!\n"
            
    if (naxis2<evtsmin or naxis2<evtscnt*imgfrac):
        is_ok=0
        msg["err"]="Too few events..aborting"
    else :
        if (gsql>0 and mod.find('both')<0) :
            ct_rate=naxis2;
            if obsexpos>0:
                ct_rate/=obsexpos
            sqlhand.write("replace into xmm_anal (obsid,mode,expos,evtraw,evtfilt,rate,xpos,ypos,done) values "+
                "(\""+obs_number+"\",\""+mod+"\",%i,%i,%i,%f,"+bore_sky_pos[mod]+",\""+today+"\");\n"%(obsexpos,evtscnt,naxis2,ct_rate))
            is_prev_sql=1
        evtscnt=naxis2

    if not os.path.exists(env["SAS_CCF"]) or os.path.isdir(env["SAS_CCF"]) or (not is_cif) or regime.find("Timing")>=0:
        return
    print "current CCF file ",env["SAS_CCF"],"\n";

    if (not os.path.exists(work_dir+"/"+source+"_"+mod+"_pat.gif") or gproc>1) :
        comp="epatplot set="+source+"_"+mod+".evt";
        #$comp+=" plotfile=".$source."_".$mod."_pat.ps";
        comp+=" plotfile="+source+"_"+mod+"_pat.gif device='/GIF'";
        comp+=" withqdp='y';"#.$source."_".$mod."_pat.qdp";
        print_sys(comp,0)

    if (gtrans>1) :
        comp="scp "+work_dir+"/"+source+"_"+mod+"_pat.gif "
        comp+=outaddr+"/pict/pile/"+obs_number+"_"+mod+".gif"
        print_sys(comp)

def gen_spectrum(tabname,colname,extexpr,specname=None,gtifile=None):
    if specname==None: specname=source+"_"+mod
    if gtifile==None: gtifile=specname
    gtifile+=".gti"
    colname_pi='TLMAX9'
    if regime.find("Timing")>=0: colname_pi='TLMAX5'
    num_channels=20479
    (max_channels,)=get_keys(tabname,[colname_pi])
    if max_channels>0 and max_channels<num_channels:
        num_channels=max_channels
#	if ($proj_channels>0 and $proj_channels<$num_channels) {$num_channels=$proj_channels;}
#	print_head "Generating $mod spectrum";
    if gproc<2 and os.path.exists(work_dir+"/"+specname+".pi"):
        print "Spectrum "+specname+".pi already exists!\n"
        return

    #my $comp="cd $work_dir/;";
    if (not os.path.exists(work_dir+"/"+specname+".pi") or gproc>1) :
        comp="evselect table="+tabname+" energycolumn="+colname+" spectrumset="+specname+".pi ";
        comp+="specchannelmin=0 specchannelmax="+num_channels+" spectralbinsize="+specbinsize+" ";
        comp+="withspecranges=yes withspectrumset=yes ";
        comp+="expression=\"(FLAG==0)&&(PATTERN<=4)";
        if os.path.exists(work_dir+"/"+gtifile):
            comp+="&&gti("+gtifile+",TIME)" 
        if extexpr:
            comp+="&&$extexpr" 
        comp+="\""
        print_sys(comp)
    
    if (not os.path.exists(work_dir+"/"+specname+".pi")):
        is_ok=0;
        msg["err"]="spectrum file "+work_dir+"/"+specname+".pi not found"

def gen_lc(tabname,lcuname,extexpr):

    lcufile=lcu_dir+"/"+lcuname+".lcu"
    comp="evselect table=$tabname rateset="+lcufile+" ";
    comp+="timebinsize=$timebinsize makeratecolumn=yes maketimecolumn=yes withrateset=yes ";
    comp+="expression=\"(FLAG==0)"
    
    if (not extexpr) or extexpr.find("PATTERN")<0: comp+="&&(PATTERN<=4)" 
    if extexpr: comp+="&&"+extexpr
    comp+="\""
    print_sys(comp)
    lc_size=0
    if os.path.exists(lcufile): (lc_size)=get_keys(lcufile,"NAXIS2",1)
    else: lc_size=-1
    if (lc_size<=0):
        is_ok=0
        if (not lc_size): os.unlink(lcufile)
        msg["err"]="Problem with lightcurve "+lcufile+" from "+tabname+" using "+extexpr+": obtained "+lc_size+"\n";
    else:
        if (gerr>2):
            sys.stderr.write("Generating lc "+lcuname+" from "+tabname+" using "+extexpr+" ["+lc_size+" lines]\n")

def gen_gti(name=None,rate=None):

    if name==None: name=source+"_"+mod+"_bkg"
    if rate==None: rate=lim_rate[mod]

#	$comp="cd $work_dir/;";
    if (not os.path.exists(work_dir+"/"+name+".gti")) or (gproc>1):
        comp="tabgtigen table="+lcu_dir+"/"+name+".lcu gtiset="+work_dir+"/"+name+".gti:STDGTI expression=\"RATE.LE.%f\""%(rate)
        print_sys(comp)
    is_ok=1
    return is_ok

##---------------------------------------------------------------------------------PREPSAT!!!
def make_gti(name,tstep=0.01):
    print_head("Removing particle background")
    if os.path.exists(lcu_dir+"/"+name+".lcu") and gproc<2:
        print "File ",lcu_dir+"/"+name+".lcu already created\n"
    else :
        gen_lc(event_list,name,"(PATTERN==0)&&(PI>10000)")

    comp="";
    #if is_ok and not os.path.exists(work_dir+"/"+name+".hst") :
        #if (exe==env["LHEASOFT"]+"/bin/fhisto") and os.path.exists(exe) :
            #comp=exe+" infile="+lcu_dir+"/"+name+".lcu outfile="+work_dir+"/"+name+".hst column=RATE binsz=0.1 clobber=yes"
        #else :
            #exe=env["HEADAS"]+"/bin/ftcopy"
            #comp=exe+" \"$lcu_dir/$name.lcu[RATE][bin RATE=::0.01 ]\" $work_dir/$name.hst clobber=yes"

        #print_sys $comp,0;
    #}

    infile=lcu_dir+"/"+name+".lcu"
    rate=pyfits.open(infile).data["RATE"]
    #irate=(rate/tstep).astype(int)
    #ipos=where(irate[1:]-irate[:-1])[0]
    from numpy import histogram,arange 
    hpos=arange(rate.min(),rate.max(),0.01)
    hdist,hpos=histogram(rate,hpos)
    
    if (is_ok): 
        gen_gti(name)
    gti_frac=1
    if os.path.exists(work_dir+"/"+name.gti) :
        rep=get_sys("ftstat "+work_dir+"/"+name+".gti+1")
        tmin,tgood=0,0
        for line in rep.split("\n"): 
            if (re.search("(S[A-Z]+)",line)):
                tgood=-tgood
            elif line.find("sum")>=0:
                mx=re.search("sum:\s+([\d\.]+)",line)
                if mx: tgood+=float(mx.group(1))
            elif line.find("min")>=0:
                mx=re.search("min:\s+([\d\.]+)",line)
                if mx and tmin==0: tmin=float(mx.group(1))
            elif line.find("max")>=0:
                mx=re.search("max:\s+([\d\.]+)",line)
                if mx: tmax=float(mx.group(1))
        if (tgood>0 and tmax>tmin) :
            gti_frac=round(tgood/float(tmax-tmin)*100.,1)
            print "GTIs cover "+gti_frac+"% of total time"
        else :
            print "GTI ERROR: %f - %f - %f"%(tmin,tmax,tgood)

def scan_level():
    global work_dir,evtscnt
    cond="obsid=\"obs_number\" and mode=\"mod\"";
    work_dir=base_dir+"/spectra"
    specname=source+"_"+mod
    specfile=work_dir+"/"+specname+"_r.pi"
    if (specfile) :
        (num_channels)=get_keys(specfile,["NAXIS2"],1)
        print "Found spectrum with "+num_channels+"bins\n"
    
    bkpset=source+"_"+mod+".rmf";
    matset=work_dir+"/"+bkpset;
    if os.path.exists(matset) :
        print "Response matrix exists\n"

    filelst=glob(base_dir+"/"+mod+"/*Evts.ds")
    if len(filelst)>0:
        (num_evts)=get_keys(filelst[0],["NAXIS2"],1)
        print "Found processed event file "+filelst[0]+" with "+num_evts+" evts\n"
        evtscnt=num_evts
        if gsql>0:
            sqlhand.write("update xmm_anal set evtraw="+evtscnt+" where "+cond+";\n")
    
    filelst=glob(work_dir+"/"+mod+".evt")
    if len(filelst)>0:
        (num_evts)=get_keys(filelst[0],["NAXIS2"],1)
        print "Found selected event file with $num_evts evts\n";
        evtscnt=num_evts
        if (gsql>0):
            sqlhand.write("update xmm_anal set evtfilt="+evtscnt+" where "+cond+";\n")


def make_square(rad,cen_x,cen_y,siz_x=0,siz_y=0,rep=0):
    sqrfrac=1;
    if (gloud>3):
        print "....making square around %i,%i [rad %f size %f,%f]\n"(cen_x,cen_y,rad,siz_x,siz_y)
    bound=0
    if (cen_x>rad+1):
        imin_x=cen_x-rad
    else:
        imin_x=1
        bound+=1
    if (cen_y>rad+1):
        imin_y=cen_y-rad
    else:
        imin_y=1
        bound+=1
    if (cen_x+rad>siz_x and siz_x>0):
        imax_x=siz_x
        bound+=1
    else:
        imax_x=cen_x+rad;
    if (cen_y+rad>siz_y and siz_y>0):
        imax_y=siz_y
        bound+=1 
    else:
        imax_y=cen_y+rad
    imin_x,imin_y,imax_x,imax_y=map(int,[imin_x,imin_y,imax_x,imax_y])
    
    if (bound==4) :
        return "" # no cut needed - small picture
    if rep==0: square=[imin_x,imax_x,imin_y,imax_y]
    else:square="[imin_x:imax_x,imin_y:imax_y]"
    sqrfrac=imax_x-imin_x
    sqrfrac*=imax_y-imin_y
    if siz_x*siz_y>0: sqrfrac/=siz_x*siz_y 
    if (sqrfrac<imgfrac): imgfrac=sqrfrac
    if (gloud>2):
        print("...selected %4.1f%% of original region \n"%(sqrfrac*100))
    if (gloud>3):
        print("....output "+square)
    return square

def get_centroid(img,ipx):
    img=img[ipx[0]:ipx[1],ipx[2]:ipx[3]]
    profx=img.mean(0)
    profx/=profx.sum()
    profy=img.mean(1)
    profy/=profy.sum()
    max_x=(prof_x*r_[ipx[0]:ipx[1]]).sum()
    max_y=(prof_y*r_[ipx[2]:ipx[3]]).sum()
    return max_x,max_y


def find_region(tabname): #
    imgfile=source+"_"+mod+".img";
    workfile=source+"_out.img";
    
    print_head("Registering extraction region")
    
    m1=re.search("(\d+),(\d+)",bore_sky_pos[mod])
    if (m1) :
        omax_x=m1.group(1)/binsize
        omax_y=m1.group(2)/binsize
    os.chdir(work_dir)
    
    if os.path.exists(imgfile) and (gproc<2):
        print "File "+imgfile+" already created\n"
        workfile=imgfile
    elif os.path.exists(workfile) and (gproc<2):#(-e $workfile and $gproc<2):
        print "Workfile "+workfile+" already created\n"
    else:
        comp="cd "+work_dir+"/;"
        comp+="evselect table="+tabname+" ";
        comp+="withimageset=yes imageset=\""+source+"_out.img\" imagebinning=\"binSize\" ";
        comp+="ximagebinsize="+binsize+" xcolumn=X ";
        comp+="yimagebinsize="+binsize+" ycolumn=Y "
        ocomp=comp
        comp+="expression=\"(FLAG==0)&&(PATTERN<=4)";
        #$comp+="&&((X,Y) IN circle($bore_sky_pos{$mod},$bore_sky_rad))";
        comp+="\";";
        print_sys(comp)
        if not os.path.exists(workfile):
            is_ok=0
            msg["err"]="Failed to extract image "+imgfile
            return -1
        
    
    (ref1,ref2,pos1,pos2,abs1,abs2,siz1,siz2)=get_keys(workfile,["CRVAL1","CRVAL2","CRPIX1","CRPIX2","CRPIX1L","CRPIX2L","NAXIS1","NAXIS2"],0)
    if gloud>2: print("read %i,%i POS %i,%i ABS %i,%i -%s- %i,%i"%(ref1,ref2,pos1,pos2,abs1,abs2,mod,omax_x,omax_y))
     
    if (ra_obj and os.path.exists(workfile)): # calculating source position on the image
        try:
            from astropy import wcs
            w=wcs.WCS(workfile)
            max_x,max_y=w.all_world2pix(ra_obj,dec_obj)
        except:
            values=print_sys("sky2xy xsky="+ra_obj+" ysky="+dec_obj+" "+workfile)
            #($values=~/Output .+\s([\d\.\+\-E]+),\s+([\d\.\+\-E]+)/s) {
    
        max_x-=abs1
        max_y-=abs2
        dist=(max_x-omax_x)**2+(max_y-omax_y)**2;
        dist= (dist>0) and sqrt(dist) or 0
        print("source [RA %4.2f ,Dec %4.2f] expected at %4.0f : %4.0f (dist %4.0f  pix)"%(ra_obj,dec_obj,max_x,max_y,dist))
        if dist<20: omax_x=max_x;omax_y=max_y; # if dist<20?
    
    dist=0
    if gproc>5: # not working at the moment
        comp="eregionanalyse srcexp='(RA,DEC) in circle(%.5f,%.5f,0.0166)'"%(ra_obj,dec_obj)
        comp+=" -V 4 imageset="+workfile
        print(comp)
    else :
        # ../3. - using smaller region for centroid calculation
        square=make_square(bore_sky_rad/binsize,max_x+abs1,max_y+abs2,siz1,siz2)
        orig_square=square
        print("...source extraction "+square)
        #better calculation using centroid
        #    $comp=$ENV{HEADAS}."/bin/ftstat \"$workfile"."$square\" centroid=yes clip=no";
        try:
            img=pyfits.open(workfile)[0].data
            ipx=square#map(int,square.replace('[','').replace(']','').split(':'))
            max_x,max_y=get_centroid(img,ipx)
            #if ($values=~/cntrd\[pixel\]: +\(([\d\.]+),([\d\.]+)\)/s) {$max_x=$1;$max_y=$2;} 
            max_x+=imin_x-1
            max_y+=imin_y-1#if ($gloud>3) {print $values."\n";}
            # here it should be rather position of maximum than the centroid?
            #if ($values=~/maximum value: +(\d+).+pixel coord: +\((\d+),(\d+)\)/s) {
            #max_v=max(img)
            top_x,top_y=unravel_index(img.argmax(), img.shape)
            max_v=img[top_x,top_y]
            if (max_v<5) :
                print "NO peak!"
                return 0
                os.unlink(workfile)
            top_x+=imin_x-1
            top_y+=imin_y-1
            dist=((max_x-top_x)**2+(max_y-top_y)**2)
            dist= (dist>0) and round(sqrt(dist),1) or 0

            print("..peak at %f,%f - barycentre distance %f"%(top_x,top_y,dist))
            if (dist>10) :
                ipx=make_square(bore_sky_rad/binsize,top_x,top_y,siz1,siz2)
                #ipx=map(int,square.replace('[','').replace(']','').split(':'))
                max_x,max_y=get_centroid(img,ipx)
                max_x+=imin_x-1
                max_y+=imin_y-1
                dist=((max_x-top_x)**2+(max_y-top_y)**2)
                dist= (dist>0) and round(sqrt(dist),1) or 0
                print("recalculated barycentre distance "+dist)
                dist=0
            avg_v=img[ipx[0]:ipx[1],ipx[2]:ipx[3]].mean()
            std_v=img[ipx[0]:ipx[1],ipx[2]:ipx[3]].std()
            if ((max_v-avg_v)<2.*std_v) :
                print "..peak too low: %5.1f  vs sigma %5.1f\n"%((max_v-avg_v),std_v)
            else :
                if (gloud>2): print "...found peak %5.1f  vs sigma %5.1f"%((max_v-avg_v),std_v)
        except:
            comp=env['LHEASOFT']+"/bin/fstatistic threshlo=INDEF threshup=INDEF $workfile";
        
        max_x-=abs1
        max_y-=abs2
        dist=(max_x-omax_x)**2+(max_y-omax_y)**2;
        if (gloud>0): print(".source found at %4.0f : %4.0f (%4.0f:%4.0f dist %4.0f pix)"%(max_x,max_y,omax_x,omax_y,dist))
    
    new_bore_sky="%.0f,%.0f"%(max_x*binsize,max_y*binsize)
    if gerr>0: 
        print "Corrected position ("+mod+") at "+new_bore_sky+" .. ("+bore_sky_rad+")"
    if (dist>max_dist_lim):
        print "..new maximum too far : %5.3f\n"%dist
        square=orig_square
    else:
        bore_sky_pos[mod]=new_bore_sky
        square=make_square(bore_sky_rad/binsize,max_x+abs1,max_y+abs2,siz1,siz2)
        #make_square($bore_sky_rad/$binsize*2,$max_x,$max_y,$siz1,$siz2);

    imgall=pyfits.open(workfile)
    imgall[0].data=imgall[0].data[ipx[0]:ipx[1],ipx[2]:ipx[3]]
    head=imgall[0].header
    head['CRPIX1']-=ipx[0]
    head['CRPIX2']-=ipx[2]
    head['NAXIS1']=ipx[1]-ipx[0]
    head['NAXIS2']=ipx[3]-ipx[2]
    imgall.writeto(imgfile)
    
    #~ exe=env['HEADAS']+"/bin/ftcopy"
    #~ if os.path.exists(exe) :
        #~ $comp="$exe \"".$workfile.$square."\" ".$imgfile;
        #~ $comp+=" clobber=yes copyall=yes";
    #~ } else {print "cannot use $exe\n"}
    #~ if ($comp and not -R $imgfile) {
    #~ # and ($workfile!=$imgfile)) {
        #~ print "..extracting : $square\n";
        #~ print_sys $comp;
        #~ if (-R $imgfile) {
            #~ unlink $workfile;
        #~ } else {
            #~ $is_ok=0;
            #~ $msg["err"]="Failed to extract subimage $imgfile";
            #~ return -1;
        #~ }
    #~ }
    
    if (gpict>2) :
        try:
            import PIL
            imax=img.max()*over
            if under: 
                imin=img.min()*under
                img-=imin
                imax-=imin
            na=img.astype("float")/imax*256
            img=imgall[0].data
            imgpng=PIL.Image.new('RGB', img.shape[::-1])
            for i in range(img.shape[0]):
              for j in range(img.shape[1]):
                imgpng.putpixel([j,img.shape[0]-i-1],tuple([int(na[i,j])]*3))
            imgpng.save(pngfile)
        except:
            print "Using to8bit conversion"
            comp=srchome+"/to8bit.bin "+imgfile+" -o tmp.fits;"
            comp+="convert -gamma 5 -sample 80x80 tmp.fits "+pngfile
            print_sys(comp)
    return 1
    
submode={}
sfilter={}

def get_bkg(event_file):
    # finds corresponding backgroud file and generates (or links) overall lightcurve .._bs_bkg.lcu
    event_list=event_file+":EVENTS"

    print_head("Background spectrum ("+mod+") from ".lc(regime))
    if (ifind(regime,"Timing")):
        gen_spectrum(event_list,"PI","(!(RAWX in ["+bore_int+"]))",source+"_"+mod+"_bkg");
    elif (ifind(regime,"Imaging") or ifind(regime,"Final")):
        (submode[mod],sfilter[mod]) = get_keys(event_file+"+1","filter","submode")
        print "extracted "+submode[mod]+sfilter[mod]+"\n"
        sub_let=(submode[mod].find("Ext")>=0) and "e" or "f"
        fil_let=(sfilter[mod].find("Med")>=0) and "m" or "t"
        mod_let=ifind(mod,"pn") and 'PN' or 'M1'
        if (mod_let=='M1' and sub_let=='e'):
            sub_let='f' # we dont have extended background from MOS detectors
        bkg_file="E1_0000"+sub_let+fil_let+"_"+mod_let+".fits"
        bkg_list=base_dir+"/"+bkg_file+":EVENTS"
        m1=re.search("\d{10}\/(.+fits)",bkg_list)
        if (is_inp and m1):
            print "Background event file used: "+m1.group(1) 
        if gerr:
            sys.stderr.write("Setting background in mode "+submode[mod]+" with filter "+sfilter[mod]+"\n")

        if os.path.exists(base_dir+"/"+bkg_file) and gproc<3:
            print "Background file ",base_dir+"/"+bkg_file," already exists\n"
        else :
            os.symlink(bkg_directory+"/"+bkg_file,base_dir+"/"+bkg_file)
            print "linking ",bkg_directory+"/"+bkg_file+"\n";
        
        if os.path.exists(base_dir+"/"+bkg_file) :
            comp="attcalc eventset=$bkg_list attitudelabel=\"fixed\" refpointlabel=\"user\" "
            comp+="fixedra="+ra_pnt+" fixeddec="+dec_pnt+" fixedposangle="+pa_pnt;
            if (ra_obj>0):
                comp+=" nominalra="+ra_obj+" nominaldec="+dec_obj
            print_sys(comp)
        
        lcu_file="E1_0000"+sub_let+fil_let+"_"+mod_let+"+lcu"
        if not os.path.exists(bkg_directory+"/"+lcu_file):
            gen_lc(bkg_list,source+"_"+mod+"_bs_bkg","(PATTERN==0)&&(PI>10000)")
            print_sys("mv "+lcu_dir+"/"+source+"_"+mod+"_bs_bkg+lcu "+bkg_directory+"/"+lcu_file)
        os.symlink(bkg_directory+"/"+lcu_file,lcu_dir+"/"+source+"_"+mod+"_bs_bkg.lcu")

        gen_gti(source+"_"+mod+"_bs_bkg")
        if os.path.exists(work_dir+"/"+source+"_"+mod+"_bs_bkg.gti"):
            os.symlink(work_dir+"/"+source+"_"+mod+"_bs_bkg+gti",work_dir+"/"+source+"_"+mod+"_bs+gti")
        gen_spectrum(bkg_list,"PI","((X,Y) IN circle("+bore_sky_pos[mod]+","+bore_sky_rad+"))",
            source+"_"+mod+"_bkg",source+"_"+mod+"_bs");
    
    if not os.path.exists(work_dir+"/"+source+"_"+mod+"_bkg.pi"):
        is_ok=0
        msg["err"]="background ("+work_dir+"/"+source+"_"+mod+"_bkg.pi) not generated"

    return bkg_list

################################################################################ main structure

def process(obs_list,entry_point=1,leave_point=6,modtxt='both',redshift=0):
    print "running events",obs_list
    #entry_point=glargs["entry_point"]
    global workdir
    rep=get_patt(entry_point,"(\d+)-(\d+)")
    if rep:
        entry_point,leave_point=map(float,rep)
    if leave_point<entry_point:
        leave_point=entry_point
    #my $bore_name;
    print "Running from lev. %i to %i ..\n"%(entry_point,leave_point)

    dirlst=("spectra","light_curves");
    if (entry_point<4) :
        dirlst.extend(modlst)

    for obs_number in (obs_list) :
        rep=init_obs(obs_number)
        if rep<0: continue
        
        for wdir in (dirlst) :
            if os.path.isdir(base_dir+"/"+wdir) :
                if gloud>1:
                    print "Directory "+wdir+" already exists\n"
                continue
            os.mkdir(base_dir+"/"+wdir)# or print "Cannot create $dir directory\n";
            if gloud>1:
                print "Created directory "+wdir+"\n" 
        
        if (gsql>0) :
            try:
                sqlhand=open(base_dir+"/output.sql","w")
            except:
                gsql=0
        os.chdir(base_dir)
        
        event_loc={}
        work_dir=base_dir+"/spectra"
        read_inp(source)
        # basic loop starts
        for step in range(entry_point,leave_point) :
            if (gloud>1) :
                rep="#"+(((gproc>0) and "=" or "-")*50)+"level step\n";
                print rep
                if glog:
                    loghand.write(rep)
            if (step==1) :
                get_data(obs_number)
                if gproc>0 and is_ok==0:
                    break
                gen_cif("ccf.cif",gproc);
                ingest_odf(1)
                print_setup()
                continue
            elif (step<=7) :
                #gen_cif("ccf.cif",0) if !$is_cif;
                if (step<=5):
                    gen_cif("ccf.cif",0)
                else:
                    gen_cif("ccf.cif",-1)
                ingest_odf(0)
            else :
                write_inp(source)
            numloc=0
            for mod in (modlst) :
                if gloud>2:
                    print "running "+mod+" ..\n";
                if (step==0):
                    scan_level()
                    continue
                    
                if (step<6 and not mod in event_loc) :
                    event_loc[mod]=find_data(mod,step)
                    if (step==2) :
                        if (not event_loc[mod] or gproc>2) :
                            proc_data(mod)
                    # not reprocessing existing data unless requested
                        continue
                if step<6 and not event_loc[mod]:
                    continue
                numloc+=1
                # we don't need eventlist after step 5
                event_file=event_loc[mod]
                event_list=event_file+":EVENTS";
                work_dir=base_dir+"/spectra";
                lcu_dir=base_dir+"/light_curves";
                os.chdir(work_dir)
                write_inp(source,mod)
                if (step==3) :
                    make_gti(source+"_"+mod+"_bkg")
                elif (step==4) :
                    print_head("Generating "+mod+" spectrum in "+regime)
                    if ifind(regime,"Final") :
                        reg_cond=""
                    else :
                        if ifind(regime,"Timing")  :
                            reg_cond="(RAWX in ["+bore_int+"])"
                        else :
                            if (find_region(event_list)>0) :
                                reg_cond="((X,Y) IN circle("+bore_sky_pos[mod]+","+bore_sky_rad+"))";
                            else:
                                event_loc[mod]=""
                                print "Aborting mode $mod"
                                if (gsql>0) :
                                    sqlhand.write("replace into xmm_anal (obsid,mode,expos,evtfilt,rate,done) values "+
                                        "(\"%s\",\"%s\",%i,-1,0,\"%s\");\n"%(obs_number,mod,obsexpos,today))
                                continue
                        gen_pileup(event_list,reg_cond);
                    
                    gen_spectrum(event_list,"PI",reg_cond);
                    gen_rmf_arf();
                
                elif (step==5): #spectra and lightcurves
                    bkg_list=get_bkg(event_file)
                    if not nh:
                        nh=get_nh(ra_obj,dec_obj)
                    calc_lc(event_list,bkg_list);
                    bin_spectrum(source+"_"+mod,3,25,0)
                elif (step==6): #xspec running
                    work_dir=base_dir+"/spectra";
                    print_load(source,mod)
                    get_rates(source)
                    print_xcm(source,mod)
                    run_xcm(source)
                    trans_data()
                elif (step==8) : #extra tests
                    comp_lc(event_list,"2,4,6,8,10")
            if numloc<1 and step>2: # no usable data
                break
            if (step==6 and mod!=modtxt) : #xspec running in combined mode
                print_head("Fitting both datasets")
                mod=modtxt;
                work_dir=base_dir+"/spectra"
                print_xcm(source,mod)
                run_xcm(source)
                trans_data()
            elif (step==7) : #transfer && cleaning
                clean_data(modlst)
            
            if is_ok==0 and gproc>0:
                break
            
            if (is_inp) :
                oufile.close()
                is_inp=0
        # basic loop ends

        if (gsql>0) :
            sqlhand.close()
            gsql=0
            rep=os.popen("grep -c \';\$\' "+base_dir+"/output.sql").read()
            rep=rep.chomp()
            print "Found "+rep+" lines\n"
            if (gtrans>1 and rep>0) :
                print_sys("cat "+base_dir+"/output.sql | ssh "+outip+" \"zmysql analyse\"")
        
        if glog:
            loghand.close()
        if not "err" in msg: 
            print "Finished\n\n"
        elif gproc>0: 
            print "Failed ("+msg["err"]+") !!!\n\n"
            if (step<5) :
                print_setup();


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-m","mode","store","modtxt")
    parser.add_option("-r","redsh","store","redshift",default=0,type="int")
    parser.add_option("-lp","lim_pn","store","lim_pn")
    parser.add_option("-lm","lim_mos","store","lim_mos")

    parser.add_option("-L","loud","store","loud")
    parser.add_option("-P","proc","store","proc")

    parser.add_option("-RA","rasc","store","ra_obj")
    parser.add_option("-DC","dec","store","dec_obj")
    
    parser.add_option("-ST","entry","store","entry_point")

    extra_args='''
    ($inp=~/redsh/i) {$redshift=$val;}
			elsif ($inp=~/name|source/i) {$source=$val;}
			elsif ($inp=~/entry|step/i) {$entry_point=$val;}
			elsif ($inp=~/leave/i) {$leave_point=$val;}
			elsif ($inp=~/lim_(pn|mos)/i) {$lim_rate{$1}=$val;}
			elsif ($inp=~/obsid/i) {$obs_number=$val;}
			elsif ($inp=~/ccf_dir/i) {$ccf_directory=$val;}
			elsif ($inp=~/bkg_dir/i) {$bkg_directory=$val;}
			elsif ($inp=~/^loud/i) {$gloud=$val;}
			elsif ($inp=~/^proc/i) {$gproc=$val;}
			elsif ($inp=~/err/i) {$gerr=$val;}
			elsif ($inp=~/pict/i) {$gpict=$val;}
			elsif ($inp=~/anal/i) {$ganal=$val;}
			elsif ($inp=~/clean/i) {$gclean=$val;}
			elsif ($inp=~/sql/i) {$gsql=$val;}
			elsif ($inp=~/check/i) {$gcheck=$val;}
			elsif ($inp=~/trans/i) {$gtrans=$val;}
			elsif ($inp=~/^rasc/i) {$ra_obj=$val;}
			elsif ($inp=~/^dec/i) {$dec_obj=$val;}
			elsif ($inp=~/^preset/i) {$gproc=1;$gsql=1;$gpict=3;$gtrans=0;$gloud=$val;}
'''
    (options, glargs) = parser.parse_args()
    process(options,**glargs)
