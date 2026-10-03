"""FIFO memory reservations for this project's native sampler qualification.

Reservations cover native attempts and their output checks. This is specific
to the qualification plan; it changes neither GPU settings nor other jobs.
"""
from collections import deque
from contextlib import contextmanager
from threading import Condition


class MemoryBudget:
    def __init__(self, capacity, event_sink=None):
        assert isinstance(capacity, int) and capacity > 0
        self.capacity = capacity
        self.condition = Condition()
        self.waiting = deque()
        self.active = {}
        self.seen = set()
        self.used = 0
        self.peak = 0
        self.events = []
        self.event_sink = event_sink
        self.aborted = None

    def record(self, **event):
        event['sequence'] = len(self.events)
        self.events.append(event)
        if self.event_sink is not None:
            try:
                self.event_sink(event)
            except BaseException:
                self.aborted = 'Reservation journal write failed'
                self.condition.notify_all()
                raise

    def abort(self, reason):
        """Stop admission after an unexpected attempt/controller exception.

        Already admitted attempts retain their limits. In particular, releasing
        a lease after an orphaned native process cannot admit another process.
        """
        with self.condition:
            if self.aborted is None:
                self.aborted = str(reason)
                self.record(action='abort', reason=self.aborted, reserved_total=self.used,
                            active_roles=len(self.active))
            self.condition.notify_all()

    @contextmanager
    def reserve(self, identifier, amount):
        assert isinstance(amount, int) and 0 < amount <= self.capacity
        with self.condition:
            assert identifier not in self.seen, 'A qualification role cannot acquire twice'
            self.seen.add(identifier)
            self.waiting.append(identifier)
            while not self.aborted and (self.waiting[0] != identifier or self.used + amount > self.capacity):
                self.condition.wait()
            if self.aborted:
                self.waiting.remove(identifier)
                self.condition.notify_all()
                raise RuntimeError('Qualification admission aborted: ' + self.aborted)
            assert self.waiting.popleft() == identifier
            self.active[identifier] = amount
            self.used += amount
            self.peak = max(self.peak, self.used)
            self.record(action='acquire', identifier=identifier, amount=amount,
                        reserved_total=self.used, active_roles=len(self.active))
            self.condition.notify_all()
        try:
            yield
        finally:
            with self.condition:
                assert self.active.pop(identifier) == amount
                self.used -= amount
                assert self.used >= 0
                self.record(action='release', identifier=identifier, amount=amount,
                            reserved_total=self.used, active_roles=len(self.active))
                self.condition.notify_all()
