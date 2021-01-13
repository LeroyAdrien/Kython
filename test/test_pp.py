from kython import pp

def test_DownloadingSequences():
        pp.DownloadingSequences('./data/Bacteria.list','./data/Archea.list','./data/refseq')==True

def test_ParsingSequences():
        assert pp.ParsingSequences('./refseq')

def test_KmerSignature():
        assert pp.KmerSignature('./data/testGenome.fna',3,False)

def test_DistanceMatrix():
        assert pp.DistanceMatrix('./refseq',3)
