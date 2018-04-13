from numpy import array


from UQSetting import UQSetting
from UQSettingFormatter import UQSettingFormatter


class Samplingresult:

    def __init__(self, uqsetting):
        """
        Create a Samplingresult using data of the given UQSetting. If
        a filename or a list of filenames is given, load the UQSetting
        from those files.
        """

        self.uq = None
        self.qoi = None
        self.time = None
        self.chunk_ref = {}

        self.onlyActive = False
        self.unitcube = False
        self.onlyValues = False
        self.raw = False

        if isinstance(uqsetting, str):
            self.loadfile(uqsetting)
        elif isinstance(uqsetting, list):
            for f in uqsetting:
                self.loadfile(f)
        else:
            self.uq = uqsetting

        self.tag_stats = self.uq.getTagsStats()
        self.sim_stats = self.uq.getSimulationStats()
        self.post_stats = self.uq.getPostprocessorStats()

    def loadfile(self, filename):
        m = UQSettingFormatter().deserializeFromFile(filename)
        uq = UQSetting.fromJson(m)
        if not self.uq:
            self.uq = uq
        else:
            self.uq.mergeStats(uq)

    def tagMatch(self, tag, chunk_ref, not_match=None):
        """
        Check if a tag matches an description.
        tag may be given as list of items or dict.
        return true if all keys in chuck_ref have the same value in tag.
        @return: false otherwise or if tag matches not.
        """
        ok = True
        if isinstance(tag, dict):
            tk = tag.items()
        else:
            tk = tag
        for keyval in chunk_ref.iteritems():
            if keyval not in tk:
                ok = False
        not_matched = False
        if not_match is not None:
            not_matched = self.tagMatch(tag, not_match)
        return ok and not not_matched

    def getSamplesSorted(self, chunktag=None, first=0, number=None,
                         not_match=None, **kwargs):
        """
        Returns the samples which match the preset parameters and are assigned
        to the chunk given by the chunktag. The samples will be ordered as
        they are in the chunk.

        @param chunktag: denotes the  chunk
        @param first: index of the first sample to return
        @param number: number of samples to return, None for all
        @param not_match: tags matching not_match will be excluded from the result
        @param kwargs: same as a dict for chunktag
        """
        chunk_ref = self.chunk_ref.copy()
        if isinstance(chunktag, dict):
            chunk_ref.update(chunktag)
        elif isinstance(chunktag, str):
            chunk_ref['type'] = chunktag
        elif chunktag is None:
            chunk_ref.update(kwargs)
        else:
            print "Warning: illegal chunktag", chunktag
            return []
        items = []
        if hasattr(self, 'index'):
            for tk, v in self.index.iteritems():
                if self.tagMatch(tk, chunk_ref, not_match):
                    for k in v:
                        t = self.tag_stats[k]
                        idx = t[0]['__index']
                        items.append((idx, k, self.sim_stats[k]))
        else:
            for k, v in self.tag_stats.iteritems():
                for t in v:
                    if self.tagMatch(t, chunk_ref, not_match):
                        idx = t['__index']
                        items.append((idx, k, self.sim_stats[k]))

        items.sort()
        out = []
        param = self.uq.getParameters()
        for i, k, v in items:
            s = [None] * 2
            if self.onlyActive:
                s[0] = param.extractActiveTuple(k)
            elif self.unitcube:
                s[0] = v['__unitcube_value']  # TODO
            else:
                s[0] = k
            r = self.post_stats[k]
            if self.raw or self.qoi is None:
                s[1] = r
            else:
                qoi = r[self.qoi]
                if self.time is not None:
                    qoi = qoi[self.time]
                else:
                    qoi = array(qoi)
                s[1] = qoi
            out.append(s)
        if self.onlyValues:
            out = [i[1] for i in out]
        else:
            out = [tuple(i) for i in out]

        if number:
            return out[first:first + number]
        return out[first:]

    def makeIndex(self):
        """
        Creates an index over the samples by tag. This speeds up searching for
        a short chunk in a huge dataset (might be a factor of 100) but this
        might not always be the case.

        The index stored in the index attribute is a dict of lists of sample
        points indexed by "atomic" tags (those reported by listChunks()) which
        are flattend into tuples of (key, value) tuples.
        """
        index = {}
        self.index = index
        for k, ts in self.tag_stats.iteritems():
            for t in ts:
                tk = self.__tagToTuple(t)
                if tk not in index:
                    index[tk] = [k]
                else:
                    index[tk].append(k)

    def setSampleFormat(self, onlyActive=False, unitcube=False, onlyValues=False, raw=False):
        """
        Specifies the format that getSamplesSorted returns.

        @param onlyActive: only add active parameters to the sample parameters TODO: Not implemented
        @param unitcube: use unitcube values (implies active subset only) TODO: not implemented
        @param onlyValues: do not return sample parameters
        @param raw: do not extract QoI values
        """
        self.onlyActive = onlyActive
        self.unitcube = unitcube
        self.onlyValues = onlyValues
        self.raw = raw

    def setQoI(self, qoi):
        """
        Set the selected QoI
        """
        self.qoi = qoi

    def setTag(self, **kwargs):
        """
        Pre-Select samples with the given tag properties.
        """
        self.chunk_ref = kwargs

    def listQoI(self):
        """
        Return the QoIs found in the dataset as a list
        """
        return self.uq.getAvailableQoI()

    def setTime(self, time):
        """
        Set a specific timestep value to return in getSamplesSorted. Default is None,
        in this case all timesteps will be returned as a vector.
        """
        self.time = time

    def listTime(self):
        """
        Return the timesteps in the dataset as a list.
        """
        # TODO: This is inefficient and kind of useless.
        return range(0, len(self.post_stats.values()[0][self.qoi]))

    def __tagToTuple(self, t):
        """
        convert a tag given as a dict to something suitable as a dict key.
        """
        return tuple(i for i in sorted(t.iteritems())
                     if not i[0].startswith('__'))

    def listChunks(self, chunk_ref=None, **kwargs):
        """
        Return the available chunks in the dataset as a list of chunk tags.
        Slow and should not be used.
        @param chunk_ref: return only those of the given prototype (Values
                           the same for all keys in prototype.)
        @param kwargs: alternate synty for chunk_ref
        """
        if chunk_ref is None:
            chunk_ref = kwargs
        tags = {}
        for ts in self.tag_stats.itervalues():
            for t in ts:
                # all key-value pairs except __index
                ok = True
                for key, val in chunk_ref.iteritems():
                    if key not in t or t[key] != val:
                        ok = False
                if ok:
                    k = self.__tagToTuple(t)
                    tags[k] = t
        return [{k: v for k, v in t.iteritems() if not k.startswith('__')}
                for t in tags.itervalues()]
