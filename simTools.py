


class simTools:
    def processRaster(spkTimes, spkGids):
        raster_dict = {'tone': {}}
        for event in raster_dict.keys(): raster_dict[event].update({'spkTimes': [], 'spkGids': []})
        for spkt_ind, spkt in enumerate(spkTimes):
            if (1500 <= (spkt) <= 5500):
                event = 'tone'
            raster_dict[event]['spkTimes'].append(spkTimes[spkt_ind])
            raster_dict[event]['spkGids'].append(spkGids[spkt_ind])
        return raster_dict

    def getPopSpks(spkTimes, spkGids, popGids, window):
        spk_gids = []
        spk_times = []
        for spk_gid_ind, spk_gid in enumerate(spkGids):
            if (spk_gid in popGids) and (window[0]<= (spkTimes[spk_gid_ind])<= window[1]):
                print(spk_gid)
                spk_gids.append(spkGids[spk_gid_ind])
                spk_times.append(spkTimes[spk_gid_ind])
        return spk_gids, spk_times






