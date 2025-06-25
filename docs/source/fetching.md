# Fetching Seismic Data

## Data Formats Overview

### MiniSEED
- Compact binary format containing only waveform samples
- Lacks station metadata (coordinates, response information)
- Contains:
  - Time series segments of fixed length
  - Basic timing and quality flags

### StationXML
- Comprehensive XML metadata standard (FDSN)
- Includes:
  - Station/channel coordinates
  - Instrument response data
  - Gain/scale information
  - Deployment histories

## Data Fetching with ObsPy

```python
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Initialize client for USP server
client = Client(base_url="http://sismo.iag.usp.br")

# Define parameters for 2025 M6.0 São Paulo event
params = {
    "network": "BR",
    "station": "*", 
    "location": "*",
    "channel": "HH?,BH?",  # All high-rate and broadband channels
    "starttime": UTCDateTime("2025-02-14T08:30:00"),
    "endtime": UTCDateTime("2025-02-14T09:00:00")
}

# Fetch waveforms
st = client.get_waveforms(**params)

print(st)
```


```{important}
1. **Time windows**: Request only necessary time spans
2. **Channel selection**: Use wildcards (HH?) for all components
3. **Metadata**: Always fetch corresponding StationXML
4. **Error handling**: Implement try/except for network issues
```

## Output Interpretation
A successful request returns an `obspy.Stream` object:
```
3 Trace(s) in Stream:
BR.????.00.HHE | 2025-02-14T08:30:00.000000Z - 2025-02-14T09:00:00.000000Z | 100.0 Hz, 180000 samples
BR.????.00.HHN | 2025-02-14T08:30:00.000000Z - 2025-02-14T09:00:00.000000Z | 100.0 Hz, 180000 samples
BR.????.00.HHZ | 2025-02-14T08:30:00.000000Z - 2025-02-14T09:00:00.000000Z | 100.0 Hz, 180000 samples
```