all: test_cgm test_preconditioner

test_cgm: cgm.fut cgm_test.fut
	futhark test cgm_test.fut

test_preconditioner: preconditioner.fut preconditioner_test.fut
	futhark test preconditioner_test.fut

.PHONY: clean

clean:
	$(RM) *.c *.actual *.expected preconditioner_test cgm_test
