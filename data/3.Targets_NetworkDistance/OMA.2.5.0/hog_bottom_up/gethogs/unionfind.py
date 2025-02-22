import collections

"""UnionFind.py

Union-find data structure. Based on Josiah Carlson's code,
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
with significant additional changes by D. Eppstein and
Adrian Altenhoff.
"""


class UnionFind(object):
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self, elements=None):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}
        if not elements is None:
            for elem in iter(elements):
                self.parents[elem] = elem
                self.weights[elem] = 1

    def __getitem__(self, obj):
        return self.find(obj)

    def find(self,obj):
        """Find and return the name of the set containing the obj."""

        # check for previously unknown obj. If unknown, add it 
        # as a new cluster
        if obj not in self.parents:
            self.parents[obj] = obj
            self.weights[obj] = 1
            return obj

        # find path of objects leading to the root
        path = [obj]
        root = self.parents[obj]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r],r) for r in roots], key=lambda x:x[0])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

    def get_components(self):
        """return a list of sets corresponding to the connected
        components of the structure."""
        comp_dict = collections.defaultdict(set)
        for elem in iter(self):
            comp_dict[self[elem]].add(elem)
        comp = list(comp_dict.values())
        return comp


