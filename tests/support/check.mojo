"""Checks that fail the test process, with a visible count."""
from std.sys import exit

struct Checker(Movable):
    var name: String
    var checks: Int
    var failures: Int

    def __init__(out self, name: String):
        self.name = name
        self.checks = 0
        self.failures = 0

    def check(mut self, condition: Bool, what: String):
        self.checks += 1
        if not condition:
            self.failures += 1
            print(self.name, "FAILED:", what)

    def finish(self):
        print(self.name, ":", self.checks, "checks,", self.failures, "failures")
        if self.failures:
            exit(1)
