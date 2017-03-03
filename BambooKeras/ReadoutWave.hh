//
// Created by herbertqiao on 3/2/17.
//

#ifndef ROOTSCRIPT_READOUTWAVE_H
#define ROOTSCRIPT_READOUTWAVE_H

#include <map>
#include <vector>
#include <TObject.h>

class ReadoutWave
        : public TObject
{
public:

    ReadoutWave()
    {}

    double trigger_offset = 0;
    double total_energy = 0;
    int begin = 0;
    int trigger = 0;
    int end = 0;
    std::map<int, std::vector<int>> detectors;

ClassDef(ReadoutWave, 1)

//    friend ostream &operator<<(ostream &out, const ReadoutWave &s)
//    {
//        out.write((char *) &s.trigger_offset, sizeof(double));
//        int channels = s.detectors.size();
//        out.write((char *) &channels, sizeof(int));
//        for (const auto &v:s.detectors) {
//            auto channel_id = v.first;
//            const auto &wave = v.second;
//            int wave_size = wave.size();
//            out.write((char *) &channel_id, sizeof(int));
//            out.write((char *) &wave_size, sizeof(int));
//            for (int j = 0; j < wave_size; j++)
//                out.write((char *) &wave[j], sizeof(int));
//        }
//        return out;
//    }
//
//    friend istream &operator>>(istream &in, ReadoutWave &s)
//    {
//        in >> s.trigger_offset;
//        int channels;
//        in >> channels;
//        for (int i = 0; i < channels; i++) {
//            int channel_id;
//            in >> channel_id;
//            int wave_size;
//            in >> wave_size;
//            vector<int> wave;
//            wave.resize(wave_size);
//            for (int j = 0; j < wave_size; i++)
//                in >> wave[j];
//            s.detectors[channel_id] = wave;
//        }
//    }

};

#endif //ROOTSCRIPT_READOUTWAVE_H
