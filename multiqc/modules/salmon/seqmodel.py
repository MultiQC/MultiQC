import numpy as np
import sys


class SeqModel:
    def __init__(self):
        self.obs3_ = None
        self.obs5_ = None
        self.exp3_ = None
        self.exp5_ = None
        self.dims_ = None
        self.valid_ = False
        self.context_length3 = 0
        self.context_length5 = 0

    def populate_model_(self, data_):
        import struct

        offset = 0
        int_struct = struct.Struct('@i')
        long_struct = struct.Struct('@q')
        double_struct = struct.Struct('@d')
        context_length = int_struct.unpack_from(data_[offset:])[0]
        offset += (int_struct.size * 3)
        offset += (int_struct.size * context_length * 3)
        nrow = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size
        ncol = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size
        offset += (nrow * ncol * double_struct.size)
        offset += (int_struct.size * 2)
        model_struct = struct.Struct('@' + 4 * context_length * 'd')
        model = model_struct.unpack_from(data_[offset:])
        model = np.array(model)
        model = model.reshape(4, context_length)
        return model

    # dname is the root directory of salmon output
    def from_file(self, dname):
        import os
        import gzip
        obs3_name = os.path.sep.join([dname, 'aux_info', 'obs3_seq.gz'])
        exp3_name = os.path.sep.join([dname, 'aux_info', 'exp3_seq.gz'])
        obs5_name = os.path.sep.join([dname, 'aux_info', 'obs5_seq.gz'])
        exp5_name = os.path.sep.join([dname, 'aux_info', 'exp5_seq.gz'])

        obs3_dat, exp3_dat = None, None
        obs5_dat, exp5_dat = None, None
        try:
            with gzip.open(obs3_name) as obs_file:
                obs3_dat = obs_file.read()
            self.obs3_ = self.populate_model_(obs3_dat)
        except IOError:
            print("Could not open file {}".format(obs3_name))
            return False

        try:
            with gzip.open(obs5_name) as obs_file:
                obs5_dat = obs_file.read()
            self.obs5_ = self.populate_model_(obs5_dat)
        except IOError:
            print("Could not open file {}".format(obs5_name))
            return False

        try:
            with gzip.open(exp3_name) as exp_file:
                exp3_dat = exp_file.read()
            self.exp3_ = self.populate_model_(exp3_dat)
        except IOError:
            print("Could not open file {}".format(exp3_name))
            return False

        try:
            with gzip.open(exp5_name) as exp_file:
                exp5_dat = exp_file.read()
            self.exp5_ = self.populate_model_(exp5_dat)
        except IOError:
            print("Could not open file {}".format(exp5_name))
            return False
        self.valid_ = True
        return True
