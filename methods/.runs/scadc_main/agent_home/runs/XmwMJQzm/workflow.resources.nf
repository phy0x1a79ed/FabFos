process {
	withName: 'p01__significance' {
		cpus = 2
		memory = { (2**(task.attempt-1)) * ('16.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('2hours' as Duration) }
	}
}
