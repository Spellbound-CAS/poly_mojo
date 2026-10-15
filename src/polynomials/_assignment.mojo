"""Call-scoped assignment and power indices shared by evaluation and substitution."""
from std.collections import Dict
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum


struct _AssignmentPlan(Movable):
    var assignments: Dict[Int64, Int]
    var powers: Dict[Tuple[Int64, UInt32], Int]
    var small_assignments: List[Tuple[Int64, Int]]
    var small_powers: List[Tuple[Int64, UInt32]]

    def __init__(out self):
        self.assignments = Dict[Int64, Int]()
        self.powers = Dict[Tuple[Int64, UInt32], Int]()
        self.small_assignments = List[Tuple[Int64, Int]]()
        self.small_powers = List[Tuple[Int64, UInt32]]()

    @always_inline
    def reserve_small(mut self, assignments: Int, powers: Int, mut budget: Budget) raises PolyError:
        var assignment_count = min(assignments, 8)
        var power_count = min(powers, 16)
        _allocation[Tuple[Int64, Int]](assignment_count, budget)
        _allocation[Tuple[Int64, UInt32]](power_count, budget)
        self.small_assignments.reserve(assignment_count)
        self.small_powers.reserve(power_count)

    @always_inline
    def assign(mut self, variable: Int64, index: Int, mut budget: Budget) raises PolyError:
        budget.charge(1)
        if not len(self.assignments):
            for entry in self.small_assignments:
                budget.charge(1)
                if entry[0] == variable:
                    raise PolyError(PolyErrorKind.DOMAIN, "duplicate assignment variable")
            if len(self.small_assignments) < 8:
                _allocation[Tuple[Int64, Int]](1, budget)
                self.small_assignments.append((variable, index))
                return
            for entry in self.small_assignments:
                _allocation[Tuple[Int64, Int]](1, budget)
                self.assignments[entry[0]] = entry[1]
        elif variable in self.assignments:
            raise PolyError(PolyErrorKind.DOMAIN, "duplicate assignment variable")
        _ = _sum(len(self.assignments), 1)
        _allocation[Tuple[Int64, Int]](1, budget)
        self.assignments[variable] = index

    @always_inline
    def assignment(self, variable: Int64, mut budget: Budget) raises PolyError -> Int:
        budget.charge(1)
        if len(self.assignments):
            return self.assignments.get(variable, -1)
        for entry in self.small_assignments:
            budget.charge(1)
            if entry[0] == variable:
                return entry[1]
        return -1

    @always_inline
    def power(mut self, pair: Tuple[Int64, UInt32], mut budget: Budget) raises PolyError -> Int:
        budget.charge(1)
        if not len(self.powers):
            for i in range(len(self.small_powers)):
                budget.charge(1)
                if self.small_powers[i][0] == pair[0] and self.small_powers[i][1] == pair[1]:
                    return i
            if len(self.small_powers) < 16:
                _allocation[Tuple[Int64, UInt32]](1, budget)
                var index = len(self.small_powers)
                self.small_powers.append(pair)
                return index
            for i in range(len(self.small_powers)):
                _allocation[Tuple[Tuple[Int64, UInt32], Int]](1, budget)
                self.powers[self.small_powers[i]] = i
        else:
            var previous = self.powers.get(pair, -1)
            if previous >= 0:
                return previous
        var index = len(self.powers)
        _ = _sum(index, 1)
        _allocation[Tuple[Tuple[Int64, UInt32], Int]](1, budget)
        self.powers[pair] = index
        # Dictionaries are lookup-only. Promotion retains first-use slot order.
        return index
