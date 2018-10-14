import numpy as np

elements_list = []


class Element:
    """

    """
    def __init__(self, el_id):
        """

        """

        self.element_id = el_id
        self.nodes = []

        elements_list.append(self)

    def add_nodes(self, node_ids):
        self.nodes = node_ids

    def ke(self, he):
        ke = np.array([1, -1], [-1, 1])
        ke *= (1/he)
        return ke
