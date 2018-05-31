<Query Kind="Program">
  <Namespace>System</Namespace>
  <Namespace>System.Collections</Namespace>
  <Namespace>System.ComponentModel</Namespace>
  <Namespace>System.Data</Namespace>
  <Namespace>System.Linq</Namespace>
  <Namespace>System.Globalization</Namespace>
</Query>

// -------------------------------------------------------------------------------------------------------

public class BatchOut
{
	public DateTime jobstart;
	public string model_file;
	public double world_size;
	public string matrix_size;
	public int seed;
	public double prep_time;
	public double end_time;
	public int maintenance;

	// sim specific param
	public double C;
	public int N;
	public double Nc;
}

public enum EndCode { DNF, FAIL, ABORT, ENDED, };

public class BatchErr
{
	public UInt64 prepsteps, simsteps;
	public double preptime, simtime;
	public EndCode endCode;

	public int RnConvErrors;
	public int PlaceProductFailed;
	public int FatalFailed;
}

public class BatchCopyNum
{
	public string[] species;
	public List<double>[] particles;
	public double interval;
}

public class BatchJob
{
	public int cores;
	public int memory;
	public List<string> run_keys;
	public String node_name;
	public String cpu_model;
	public DateTime jobstart;
	public DateTime jobend;
	public int RunsTerminated;
	public int RunsFailed;
}

// -------------------------------------------------------------------------------------------------------

Dictionary<string, BatchOut> results1 = new Dictionary<string, BatchOut>();
Dictionary<string, BatchErr> results2 = new Dictionary<string, BatchErr>();
Dictionary<string, BatchJob> results3 = new Dictionary<string, BatchJob>();
Dictionary<string, BatchCopyNum> results4 = new Dictionary<string, BatchCopyNum>();

// -------------------------------------------------------------------------------------------------------

public string GetKey(string file)
{
	var name = Path.GetFileNameWithoutExtension(file);
	int pos1 = name.IndexOf('_');
	return name.Substring(pos1 + 1);    //.PadLeft(2,'_');
}

void ParseOutputFile(string file)
{
	var key = GetKey(file);
	var o = new BatchOut();
	foreach (var line in File.ReadAllLines(file))
	{
		var s = line.Split(new[] { "=>", "=", " [" }, StringSplitOptions.RemoveEmptyEntries);
		if (s.Length == 0) { continue; }
		var lt = s[0].Trim();

		if (lt == "time local") { var sub = line.Substring(13); o.jobstart = DateTime.ParseExact(sub, "ddd MMM d HH:mm:ss yyyy", new CultureInfo("en-US"), DateTimeStyles.AllowWhiteSpaces); }
		if (lt == "model file") o.model_file = s[1];

		if (lt == "world size") o.world_size = Double.Parse(s[1]);
		if (lt == "matrix size") o.matrix_size = s[1].Trim();
		if (lt == "seed") o.seed = Int32.Parse(s[1].Substring(3), NumberStyles.HexNumber);
		if (lt == "prep_time") o.prep_time = Double.Parse(s[1]);
		if (lt == "end_time") o.end_time = Double.Parse(s[1]);
		if (lt == "maintenance") o.maintenance = Int32.Parse(s[1]);

		if (s.Length == 3) // parse sim specific variables here
		{
			if (lt == "C") o.C = Double.Parse(s[2]);
			if (lt == "N") o.N = Int32.Parse(s[2]);
			if (lt == "Nc_th") o.Nc = Double.Parse(s[2]);
		}
	}
	results1[key] = o;
}

void ParseErrorFile(string file)
{
	var key = GetKey(file);
	var e = new BatchErr();
	foreach (var line in File.ReadAllLines(file))
	{
		var s = line.Split(new[] { "] : " }, StringSplitOptions.RemoveEmptyEntries);
		if (s.Length != 2) { Console.WriteLine("unknown: " + line); continue; }
		var lt = s[1].Trim();

		if (lt.StartsWith("Start of ")) { }
		else if (lt.StartsWith("Model simulation")) { }
		else if (lt.StartsWith("Rn didn't converge")) e.RnConvErrors++;
		else if (lt.Contains("placing product(s) failed")) e.PlaceProductFailed++;
		else if (lt.Contains("placing product failed")) e.PlaceProductFailed++;
		else if (s[0].Contains("[FATAL]")) e.FatalFailed++;
		else if (lt.StartsWith("Pre-simulation "))
		{
			var split = lt.Split(new[] { ' ' }, StringSplitOptions.None);
			e.preptime = Double.Parse(split[3]);
			e.prepsteps = UInt64.Parse(split[6]);
		}
		else if (lt.StartsWith("Simulation "))
		{
			var split = lt.Split(new[] { ' ' }, StringSplitOptions.None);
			int pos = 3;
			if (split[1] == "failed") { pos++; e.endCode = EndCode.FAIL; }
			else if (split[1] == "aborted") { pos++; e.endCode = EndCode.ABORT; }
			else e.endCode = EndCode.ENDED;
			e.simtime = Double.Parse(split[pos]);
			e.simsteps = UInt64.Parse(split[pos + 3]);
		}
		else Console.WriteLine("unknown: " + line);
	}
	results2[key] = e;
}

void ParseCopyNumFile(string file)
{
	var key = GetKey(file);
	var c = new BatchCopyNum();
	foreach (var line in File.ReadAllLines(file))
	{
		var s = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
		if (c.species == null)
		{
			c.species = s.Skip(1).ToArray();
			c.particles = new List<double>[s.Length - 1];
			for (var i = 1; i < s.Length; i++) c.particles[i - 1] = new List<double>();
		}
		else
		{
			if (c.interval == 0) c.interval = Double.Parse(s[0]);
			for (var i = 1; i < s.Length; i++) c.particles[i - 1].Add(Double.Parse(s[i]));
		}
	}
	results4[key] = c;
}

void ParseJobOutputFile(string file)
{
	var key = file.Substring(file.LastIndexOf('.') + 2);
	BatchJob j;
	if (results3.ContainsKey(key)) j = results3[key]; else results3[key] = j = new BatchJob();
	j.run_keys = new List<string>();
	foreach (var line in File.ReadAllLines(file))
	{
		if (line.StartsWith("Cores")) j.cores = Int32.Parse(line.Split(new[] { ' ' })[1]);
		else if (line.StartsWith("CPU")) j.cpu_model = line.Substring(3).Trim();
		else if (line.StartsWith("Run ")) j.run_keys.Add(line.Substring(line.IndexOf("key") + 4));
		else if (line.Contains("started")) j.jobstart = DateTime.ParseExact(line.Substring(line.IndexOf(" on ") + 4), "ddd MMM d HH:mm:ss CEST yyyy", new CultureInfo("en-US"));
		else if (line.Contains("ended")) j.jobend = DateTime.ParseExact(line.Substring(line.IndexOf(" on ") + 4), "ddd MMM d HH:mm:ss CEST yyyy", new CultureInfo("en-US"));
		else if (line.Contains("Host")) j.node_name = line.Split(new[] { ' ' })[1].Trim();
		else if (line.StartsWith("Memory")) j.memory = Int32.Parse(line.Split(new[] { ' ' })[1]);
		else Console.WriteLine("unknown: " + line);
	}
}

void ParseJobErrorFile(string file)
{
	var key = file.Substring(file.LastIndexOf('.') + 2);
	BatchJob j;
	if (results3.ContainsKey(key)) j = results3[key]; else results3[key] = j = new BatchJob();
	foreach (var line in File.ReadAllLines(file))
	{
		if (line.StartsWith("stopos")) { }
		else if (line.Contains("Terminated")) j.RunsTerminated++;
		else if (line.Contains("Segmentation")) j.RunsFailed++;
		else Console.WriteLine("unknown: " + line);
	}
}

// -------------------------------------------------------------------------------------------------------

void Main()
{
	// Process
	//var dir = @"\\storage01\data\Amolf\projects\tenwolde-projects\Engineering\eGFRD\Resources\Lisa\zx84_compare_org_new\output_o3_fix\";
	var dir = @"d:\Data\eGFRD\Equil\";
	foreach (var file in Directory.EnumerateFiles(dir))
	{
		var ext = Path.GetExtension(file);
		var name = Path.GetFileNameWithoutExtension(file);
		if (ext.StartsWith(".out")) ParseOutputFile(file);
		if (ext.StartsWith(".err")) ParseErrorFile(file);
		if (ext.StartsWith(".cn")) ParseCopyNumFile(file);
		if (name == "gfrd_sweep.job" && ext.StartsWith(".o")) ParseJobOutputFile(file);
		if (name == "gfrd_sweep.job" && ext.StartsWith(".e")) ParseJobErrorFile(file);
	}

	// Write Results
	using (var sw = new StreamWriter(Path.Combine(dir, "results.csv")))
	{
		sw.Write("Key\tjobstart\tmodel_file\tworld_size\tmatrix_size\tseed\tprep_time\tend_time\tmaintenance");
		sw.Write("\tC\tN\tNc");
		sw.Write("\tPrepSteps\tPrepTime\tSimSteps\tSimTime\tEndCode\tRnConvErrors\tPlaceProductFailed\tFailed");
		sw.Write("\tNavg");
		sw.WriteLine();

		foreach (var key in results1.Keys.Concat(results2.Keys).Distinct().OrderBy(k => k))
		{
			sw.Write(key); sw.Write("\t");
			if (results1.Keys.Contains(key))
			{
				var o = results1[key];
				sw.Write(o.jobstart); sw.Write("\t");
				sw.Write(o.model_file); sw.Write("\t");
				sw.Write(o.world_size); sw.Write("\t");
				sw.Write(o.matrix_size); sw.Write("\t");
				sw.Write("0x"); sw.Write(o.seed.ToString("X8")); sw.Write("\t");
				sw.Write(o.prep_time); sw.Write("\t");
				sw.Write(o.end_time); sw.Write("\t");
				sw.Write(o.maintenance); sw.Write("\t");

				sw.Write(o.C); sw.Write("\t");
				sw.Write(o.N); sw.Write("\t");
				sw.Write(o.Nc); sw.Write("\t");
			}
			else sw.Write(new string('\t', 11));

			if (results2.Keys.Contains(key))
			{
				var e = results2[key];
				sw.Write(e.prepsteps); sw.Write("\t");
				sw.Write(e.preptime); sw.Write("\t");
				sw.Write(e.simsteps); sw.Write("\t");
				sw.Write(e.simtime); sw.Write("\t");
				sw.Write(e.endCode); sw.Write("\t");
				sw.Write(e.RnConvErrors); sw.Write("\t");
				sw.Write(e.PlaceProductFailed); sw.Write("\t");
				sw.Write(e.FatalFailed); sw.Write("\t");
			}
			else sw.Write(new string('\t', 8));

			if (results4.Keys.Contains(key))
			{
				var c = results4[key];  // parse sim specific variables here
				sw.Write(c.particles != null ? c.particles[2].Average() : 0); sw.Write("\t");
			}
			else sw.Write(new string('\t', 1));

			sw.WriteLine();
		}

		sw.WriteLine();
		sw.WriteLine();
		sw.WriteLine("Job\tcores\tmemory\tnode_name\tcpu_model\tstart\tend\tTerminated\tFailed\tkeys");
		foreach (var key in results3.Keys.OrderBy(k => k))
		{
			sw.Write(key); sw.Write("\t");
			var j = results3[key];
			sw.Write(j.cores); sw.Write("\t");
			sw.Write(j.memory); sw.Write("\t");
			sw.Write(j.node_name); sw.Write("\t");
			sw.Write(j.cpu_model); sw.Write("\t");
			sw.Write(j.jobstart); sw.Write("\t");
			sw.Write(j.jobend); sw.Write("\t");
			sw.Write(j.RunsTerminated); sw.Write("\t");
			sw.Write(j.RunsFailed); sw.Write("\t");
			if (j.run_keys!=null) foreach (var k in j.run_keys) { sw.Write(k); sw.Write("\t"); }
			sw.WriteLine();
		}
		sw.WriteLine();
		sw.WriteLine();
	}
}