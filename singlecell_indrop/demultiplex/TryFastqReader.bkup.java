package demultiplex;

import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

public class TryFastqReader {

	public static void main(String[] args) throws IOException {
		System.out.println(System.currentTimeMillis());
		
		List<String> cbs = Files.readAllLines(FileSystems.getDefault().getPath("/PHShome/si992/commonscripts/singlecell_indrop/WT1_barcodes.txt"));
		HashMap<String, ArrayList<FastqRecord>> barcodes_and_fastqs = new HashMap<String, ArrayList<FastqRecord>>();
		ArrayList<FastqRecord> invalid_fastqs = new ArrayList<FastqRecord>();
		HashMap<String, ArrayList<String>> barcodes_and_umis = new HashMap<String, ArrayList<String>>();
		
		// Initialize hashmap
		Iterator<String> it = cbs.iterator();
		while (it.hasNext())
		{
			String bc = it.next();
			barcodes_and_fastqs.put(bc, new ArrayList<FastqRecord>());
			barcodes_and_umis.put(bc, new ArrayList<String>());
		}
		
		FastqReader freader = new FastqReader(new File("/data/langenau/singlecell_prkdc/processing/demultiplex_splits2/Zebrafish-WT1.processed.fastq.gz.split/00000"));
		while (freader.hasNext())
		{
			FastqRecord frecord = freader.next();
			String readname = frecord.getReadHeader();
			String[] splits = readname.split(":");
			String bc = splits[7].replaceAll("CELL_", "");
			String umi = splits[8].replaceAll("UMI_", "");
			if (cbs.contains(bc)) {
				barcodes_and_fastqs.get(bc).add(frecord);
				barcodes_and_umis.get(bc).add(umi);
			}  else {
				invalid_fastqs.add(frecord);
			}
			
		}
		freader.close();
		
		
		Iterator<String> bc_iterator = barcodes_and_fastqs.keySet().iterator();
		while (bc_iterator.hasNext())
		{
			String thisbc = bc_iterator.next();
			FastqWriter fw = new FastqWriterFactory().newWriter(new File("/PHShome/si992/RECYCLE_BIN/out_"+thisbc+".fq"));
			FileWriter uw = new FileWriter("/PHShome/si992/RECYCLE_BIN/out_"+thisbc+".umi");
			ArrayList<FastqRecord> fqrecords = barcodes_and_fastqs.get(thisbc);
			ArrayList<String> umirecords = barcodes_and_umis.get(thisbc);
			Iterator<FastqRecord> fqiterator = fqrecords.iterator();
			Iterator<String> umiiterator = umirecords.iterator();
			while (fqiterator.hasNext())
			{
				fw.write(fqiterator.next());
			}
			while (umiiterator.hasNext())
                        {
                                uw.write(umiiterator.next() + "\n");
                        }
			fw.close();
			uw.close();
		}
		FastqWriter invalid_fw = new FastqWriterFactory().newWriter(new File("/PHShome/si992/RECYCLE_BIN/invalid_fastqs.fq"));
		Iterator<FastqRecord> invalidfqiterator = invalid_fastqs.iterator();
		while (invalidfqiterator.hasNext())
		{
			invalid_fw.write(invalidfqiterator.next());
		}
		invalid_fw.close();
		
		System.out.println(System.currentTimeMillis());
	}


}
