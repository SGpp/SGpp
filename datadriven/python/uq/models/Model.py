class Model(object):

    def __init__(self):
        self.num_results = 1

    def evaluate(self, sample):
        raise NotImplementedError

    def evaluate_set(self, samples, distributed=False):
        ans = np.ndarray(num_results, len(samples))
        for i, sample in enumerate(samples):
            ans[:, i] = self.evaluate(sample)
        return ans
