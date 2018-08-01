package NovelFitters;

public enum LundPID {
	Pion(211, 0.13957),
	Kaon(321, 0.493677),
	Proton(2212, 0.938272),
	Electron(11,0.000511);
	private final int m_lundCode;
	private final double m_mass;
	
	LundPID(int pid, double mass)
	{
		this.m_lundCode=pid;
		this.m_mass=mass;
	}
	public int lundCode()
	{
		return m_lundCode;	
	
	}
	public double mass()
	{
		return m_mass;	
	}
	
}
