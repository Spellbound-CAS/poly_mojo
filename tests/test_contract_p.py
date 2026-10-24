"""Check the frozen Contract P declarations without changing the build tooling.

Run from the repository root: pixi run python tests/test_contract_p.py
The snapshot is reviewed release data; this test has no update mode.
"""
import ast
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from mojo_api import parse  # noqa: E402

SNAPSHOT = ROOT / "tests/contracts/contract-p.json"
REVISION_ONE = ROOT / "tests/contracts/contract-p-v1.json"
NUMERICAL = ROOT / "tests/contracts/contract-numeric.json"
NUMERICAL_ONE = ROOT / "tests/contracts/contract-numeric-v1.json"
DESCARTES_STREAM = {"DescartesRootItem", "DescartesRootSink", "DescartesStreamChecker"}
DESCARTES_DEFAULT = {"real_roots_exact", "refine_real_roots_exact", "check_real_roots_exact", "real_roots_exact_to", "real_root_multiplicities", "check_real_root_multiplicities", "RootMultiplicityCertificate", "RealRootStreamChecker"}
FINAL_REVISION_TWO = ROOT / "tests/contracts/contract-p-v2.json"
INITIAL_REVISION_TWO = ROOT / "tests/contracts/contract-p-v2-initial.json"
REVISION_THREE = ROOT / "tests/contracts/contract-p-v3.json"
REVISION_FOUR = ROOT / "tests/contracts/contract-p-v4.json"
REVISION_FIVE = ROOT / "tests/contracts/contract-p-v5.json"
REVISION_SIX = ROOT / "tests/contracts/contract-p-v6.json"
REVISION_SEVEN = ROOT / "tests/contracts/contract-p-v7.json"
REVISION_EIGHT = ROOT / "tests/contracts/contract-p-v8.json"
REVISION_NINE = ROOT / "tests/contracts/contract-p-v9.json"
REVISION_TEN = ROOT / "tests/contracts/contract-p-v10.json"
REVISION_ELEVEN = ROOT / "tests/contracts/contract-p-v11.json"
REVISION_TWELVE = ROOT / "tests/contracts/contract-p-v12.json"
REVISION_THIRTEEN = ROOT / "tests/contracts/contract-p-v13.json"
REVISION_FOURTEEN = ROOT / "tests/contracts/contract-p-v14.json"
# Owner: "break callers, use new API. change callers of old API to new api".
FINAL_CAD = {"CADFinalCell", "CADFinalCertificate", "cad_final", "check_cad_final"}
PACKED = {"PackedPoly", "packed_mul", "packed_pow", "packed_divide_exact"}
DESCARTES = {"DescartesLeaf", "DescartesCertificate", "DescartesRootResult", "check_descartes", "descartes_intervals"}


def declaration(item):
    return {
        "name": item.name,
        "signature": item.signature,
        "members": [declaration(member) for member in item.members],
    }


def current_contract(*, include_numerical=False):
    exports = []
    numerical_names = {item["name"] for item in json.loads(NUMERICAL.read_text())["exports"]}
    tree = ast.parse((ROOT / "src/polynomials/__init__.mojo").read_text())
    modules = {}
    for node in tree.body:
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.level != 1 or any(alias.asname for alias in node.names):
            raise ValueError("review new root export syntax before changing Contract P")
        path = ROOT / f"src/polynomials/{node.module}.mojo"
        if node.module not in modules:
            modules[node.module] = {item.name: item for item in parse(path)[1]}
        for alias in node.names:
            if not include_numerical and alias.name in numerical_names:
                continue
            exports.append({
                "module": node.module,
                **declaration(modules[node.module][alias.name]),
            })
    return {"contract": "P", "revision": 16, "exports": sorted(exports, key=lambda x: x["name"])}


SQRT = {
    "name": "sqrt",
    "signature": "def sqrt(self, mut budget: Budget) raises PolyError -> Self",
    "members": [],
}


TRACE_FIELD = {
    "name": "trace",
    "signature": "var trace: Optional[GroebnerTrace[Self.C]]",
    "members": [],
}

TRACE_CONSTRUCTORS = [
    {
        "name": "__init__",
        "signature": "def __init__(out self, var order: MonomialOrder, var basis: List[SparsePoly[Self.C]], var cofactors: Optional[List[List[SparsePoly[Self.C]]]])",
        "members": [],
    },
    {
        "name": "__init__",
        "signature": "def __init__(out self, var order: MonomialOrder, var basis: List[SparsePoly[Self.C]], var cofactors: Optional[List[List[SparsePoly[Self.C]]]], var trace: Optional[GroebnerTrace[Self.C]])",
        "members": [],
    },
]


def field(name, kind):
    return {"name": name, "signature": f"var {name}: {kind}", "members": []}


# Revision 9 appends certificate fields for witness-based CAD proofs.
WITNESS_FIELDS = {
    "AlgebraicPolynomialEvaluation": [field("witness", "List[RealAlgebraic]"), field("representations", "List[DensePoly[QQ]]")],
    "AlgebraicRootsCertificate": [field("filters", "List[AlgebraicRootFilter]"), field("element", "RealAlgebraic"),
                                  field("representations", "List[DensePoly[QQ]]")],
    "CADCell": [field("witnesses", "List[CADSignWitness]")],
}
WITNESS = {"AlgebraicRootFilter", "CADSignWitness"}
# Revision 10: exact real values and their certificates hold Descartes
# evidence instead of Sturm sequences.
DESCARTES_VALUES = {"AlgebraicDescartesCertificate", "check_algebraic_descartes"}
DESCARTES_FIELDS = {
    ("AlgebraicComparisonCertificate", "shared"): "var shared: Optional[RootInterval]",
    ("AlgebraicRootsCertificate", "all_roots"): "var all_roots: AlgebraicDescartesCertificate",
    ("AlgebraicArithmeticCertificate", "selection"): "var selection: Optional[RootInterval]",
}
# Revision 10 adds the Gröbner basis check.
GROEBNER_CHECK = {"check_groebner_basis"}

# Revision 10 appends sign witnesses to streamed CAD terminal events.
STREAM_WITNESS_FIELDS = {"CADStreamEvent": [field("witnesses", "List[CADSignWitness]")]}

# Revision 10 (owner: "Approve the full coverings proposal") adds conflict-driven
# coverings evidence for sign conjunctions, and CADFiber methods producing and
# checking it; poly.cad.Solution.evidence() returns CADCoveringsCertificate.
COVERINGS = {"CoveringBound", "CoveringInterval", "CoveringNode", "CADCoveringsCertificate", "coverings_solve", "check_coverings"}
COVERINGS_FIBER = [
    {"name": "search", "signature": "def search(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CADCoveringsCertificate", "members": []},
    {"name": "check", "signature": "def check(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CADCoveringsCertificate, mut budget: Budget) raises PolyError -> Bool", "members": []},
]

# Revision 11 (owner: "Approve revision 11 with (a)") adds the Descartes root
# count; poly's RootCount.evidence() returns DescartesCountCertificate.
COUNT = {"DescartesCountCertificate", "descartes_count", "check_descartes_count"}

# Revision 12 (owner: "Approve sign witness kind 4") adds CADSignWitness kind
# four, a section of a dividing top-level basis member; no declaration changes.
# Revision 13 (owner: "I approve", to the reduced McCallum projection) adds only
# coefficients from the leading one down to the first nonzero constant; no
# declaration changes.

def migrated_declaration(item, *, storage=True, certify=True, sqrt=True, trace=True, witnesses=True):
    """Allow only the reviewed changes in older snapshots: the revision-3 dense
    storage change (when `storage`), the revision-4 `certify` option (when
    `certify`), the revision-5 `TruncatedSeries.sqrt` method (when `sqrt`) and
    the revision-6 trace (when `trace`): `GroebnerBasis.trace` with explicit
    constructors that keep the three-argument form, and `groebner`'s `trace`
    option. The exports `divide_exact_dense` (revision 5), `TraceTerm` and
    `GroebnerTrace` (revision 6), `PackedPoly`, `packed_mul`, `packed_pow` and
    `packed_divide_exact` (revision 7) and the Descartes certificate (revision 8)
    are new, so no older baseline holds them. Revision 9 appends witness fields
    to three CAD certificates, in every baseline before it (when `witnesses`),
    and adds the exports `AlgebraicRootFilter` and `CADSignWitness`. Revision 10
    removes `RealAlgebraic.sequence` and retypes the Sturm fields in
    DESCARTES_FIELDS, in every older baseline, appends `witnesses` to
    `CADStreamEvent` in every older baseline, and adds the exports in
    DESCARTES_VALUES, GROEBNER_CHECK and COVERINGS and the CADFiber methods
    `search` and its `check` overload."""
    item = copy.deepcopy(item)
    if item["name"] == "RealAlgebraic":
        assert any(member["name"] == "sequence" for member in item["members"])
        item["members"] = [member for member in item["members"] if member["name"] != "sequence"]
    for member in item["members"]:
        if (item["name"], member["name"]) in DESCARTES_FIELDS:
            member["signature"] = DESCARTES_FIELDS[(item["name"], member["name"])]
    if item["name"] == "CADFiber":
        members = item["members"]
        assert all(member["name"] != "search" for member in members)
        members[next(i for i, member in enumerate(members) if member["name"] == "solve") + 1:0] = copy.deepcopy(COVERINGS_FIBER)
    if item["name"] in STREAM_WITNESS_FIELDS:
        assert all(member["name"] != "witnesses" for member in item["members"])
        item["members"].extend(copy.deepcopy(STREAM_WITNESS_FIELDS[item["name"]]))
    if witnesses and item["name"] in WITNESS_FIELDS:
        names = {member["name"] for member in item["members"]}
        assert names.isdisjoint(member["name"] for member in WITNESS_FIELDS[item["name"]])
        item["members"].extend(copy.deepcopy(WITNESS_FIELDS[item["name"]]))
    if trace and item["name"] == "GroebnerBasis":
        decorator = "@fieldwise_init\n"
        assert item["signature"].startswith(decorator)
        item["signature"] = item["signature"][len(decorator):]
        members = item["members"]
        assert all(member["name"] != "trace" for member in members)
        members.insert(next(i for i, member in enumerate(members) if member["name"] == "cofactors") + 1, copy.deepcopy(TRACE_FIELD))
        members.extend(copy.deepcopy(TRACE_CONSTRUCTORS))
    if sqrt and item["name"] == "TruncatedSeries":
        members = item["members"]
        assert all(member["name"] != "sqrt" for member in members)
        members.insert(next(i for i, member in enumerate(members) if member["name"] == "exp") + 1, copy.deepcopy(SQRT))
    if certify and item["name"] == "groebner":
        old = "track_cofactors: Bool = False) raises"
        assert old in item["signature"]
        item["signature"] = item["signature"].replace(old, "track_cofactors: Bool = False, certify: Bool = True) raises")
    if trace and item["name"] == "groebner":
        old = "certify: Bool = True) raises"
        assert old in item["signature"]
        item["signature"] = item["signature"].replace(old, "certify: Bool = True, trace: Bool = False) raises")
    if storage and item["name"] == "DensePoly":
        members = item["members"]
        field = next(member for member in members if member["name"] == "coefficients")
        assert field["signature"] == "var coefficients: List[Self.C]"
        field["signature"] = "var coefficients: CoefficientBuffer[Self.C]"
        factory = next(i for i, member in enumerate(members) if member["name"] == "from_coefficients")
        added = copy.deepcopy(members[factory])
        added["signature"] = added["signature"].replace("List[Self.C]", "CoefficientBuffer[Self.C]")
        members.insert(factory + 1, added)
    return item


class ContractP(unittest.TestCase):
    def test_public_declarations_are_frozen(self):
        expected = json.loads(SNAPSHOT.read_text())
        self.maxDiff = None
        self.assertEqual(current_contract(), expected,
                         "Contract P changed: agree the change and migration with consumers first")

    def test_numerical_declarations_are_separately_frozen(self):
        expected = json.loads(NUMERICAL.read_text())
        names = {item["name"] for item in expected["exports"]}
        actual = [item for item in current_contract(include_numerical=True)["exports"]
                  if item["name"] in names]
        self.assertEqual(actual, expected["exports"])
        self.assertTrue(names.isdisjoint(item["name"] for item in json.loads(SNAPSHOT.read_text())["exports"]))

    def test_numerical_revision_one_has_only_reviewed_descartes_changes(self):
        baseline = json.loads(NUMERICAL_ONE.read_text())
        current = {item["name"]: item for item in json.loads(NUMERICAL.read_text())["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, DESCARTES_STREAM)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                if item["name"] not in DESCARTES_DEFAULT:
                    self.assertEqual(item, current[item["name"]])

    def test_final_revision_two_has_only_reviewed_storage_change(self):
        baseline = json.loads(FINAL_REVISION_TWO.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertEqual(migrated_declaration(item), current[item["name"]])

    def test_revision_fourteen_has_only_approved_final_cad_additions(self):
        baseline = json.loads(REVISION_FOURTEEN.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertEqual(item, current[item["name"]])

    def test_exports_are_unique(self):
        names = [item["name"] for item in current_contract(include_numerical=True)["exports"]]
        self.assertEqual(len(names), len(set(names)))

    def test_revision_one_has_only_reviewed_storage_change(self):
        baseline = json.loads(REVISION_ONE.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item), current[item["name"]])

    def test_initial_revision_two_has_only_reviewed_storage_change(self):
        baseline = json.loads(INITIAL_REVISION_TWO.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item), current[item["name"]])

    def test_revision_three_has_only_reviewed_certify_option(self):
        baseline = json.loads(REVISION_THREE.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False), current[item["name"]])

    def test_revision_four_has_only_reviewed_additions(self):
        baseline = json.loads(REVISION_FOUR.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, {"divide_exact_dense", "GroebnerTrace", "TraceTerm"} | PACKED | DESCARTES | WITNESS | DESCARTES_VALUES | GROEBNER_CHECK | COVERINGS | COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False, certify=False), current[item["name"]])

    def test_revision_five_has_only_reviewed_trace(self):
        baseline = json.loads(REVISION_FIVE.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, {"GroebnerTrace", "TraceTerm"} | PACKED | DESCARTES | WITNESS | DESCARTES_VALUES | GROEBNER_CHECK | COVERINGS | COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False, certify=False, sqrt=False), current[item["name"]])

    def test_revision_six_has_only_reviewed_additions(self):
        baseline = json.loads(REVISION_SIX.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, PACKED | DESCARTES | WITNESS | DESCARTES_VALUES | GROEBNER_CHECK | COVERINGS | COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False, certify=False, sqrt=False, trace=False), current[item["name"]])

    def test_revision_seven_has_only_reviewed_descartes_additions(self):
        baseline = json.loads(REVISION_SEVEN.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, DESCARTES | WITNESS | DESCARTES_VALUES | GROEBNER_CHECK | COVERINGS | COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False, certify=False, sqrt=False, trace=False), current[item["name"]])

    def test_revision_eight_has_only_reviewed_witnesses(self):
        baseline = json.loads(REVISION_EIGHT.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, WITNESS | DESCARTES_VALUES | GROEBNER_CHECK | COVERINGS | COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False, certify=False, sqrt=False, trace=False), current[item["name"]])

    def test_revision_nine_has_only_reviewed_revision_ten_changes(self):
        baseline = json.loads(REVISION_NINE.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, DESCARTES_VALUES | GROEBNER_CHECK | COVERINGS | COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertIn(item["name"], current)
                self.assertEqual(migrated_declaration(item, storage=False, certify=False, sqrt=False, trace=False, witnesses=False),
                                 current[item["name"]])

    def test_revision_ten_has_only_reviewed_count_additions(self):
        baseline = json.loads(REVISION_TEN.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - {item["name"] for item in baseline["exports"]}, COUNT | FINAL_CAD)
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertEqual(item, current[item["name"]])

    def test_revision_eleven_has_no_declaration_changes(self):
        baseline = json.loads(REVISION_ELEVEN.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - FINAL_CAD, {item["name"] for item in baseline["exports"]})
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertEqual(item, current[item["name"]])

    # Revision 14 (owner: "Approve both", to the per-case projection choice and
    # Brown's rule; this revision is the first) makes the per-case projection
    # choice cad.decompose's default; cad.decompose is outside the Contract P
    # declarations, so none changes.
    def test_revision_thirteen_has_no_declaration_changes(self):
        baseline = json.loads(REVISION_THIRTEEN.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - FINAL_CAD, {item["name"] for item in baseline["exports"]})
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertEqual(item, current[item["name"]])

    def test_revision_twelve_has_no_declaration_changes(self):
        baseline = json.loads(REVISION_TWELVE.read_text())
        current = {item["name"]: item for item in current_contract()["exports"]}
        self.assertEqual(set(current) - FINAL_CAD, {item["name"] for item in baseline["exports"]})
        for item in baseline["exports"]:
            with self.subTest(name=item["name"]):
                self.assertEqual(item, current[item["name"]])


if __name__ == "__main__":
    unittest.main()
