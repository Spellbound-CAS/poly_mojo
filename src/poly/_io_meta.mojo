"""Private static metadata for the value interchange dispatcher."""
trait _Record(Copyable, Deinitable):
    comptime _Format: Copyable & Deinitable
    comptime _Item: Copyable & Deinitable
    comptime _Power: Copyable & Deinitable

struct _Dense(Copyable, Movable): pass
struct _Sparse(Copyable, Movable): pass
struct _Series(Copyable, Movable): pass
struct _Basis(Copyable, Movable): pass
struct _Fraction(Copyable, Movable): pass
struct _Batch(Copyable, Movable): pass
