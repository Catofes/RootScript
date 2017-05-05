###Convert Detail.

####Get ReadoutWave. 

First of all, we should be able to get ReadoutWave directly from our readout plane. That should be something like: 

```
// First int means the channel_id, second is the ReadoutWave for that channel.
map<int, OneReadoutWave> readout_waves
```

In our design, each channel will record with 5MHz and the time window is 512 samples, 512x0.2us=102.4us. 
For simplification, I regard that the channel will record the electrons count which hit the channel.
So the OneReadoutWave structure should be something like:

```
//The size of the vector should be exactly 512.
typedef std::vector<int> OneReadoutWave
```