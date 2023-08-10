TEST_FRAGMENTS = \
	lib/github.com/jonesz/krylov-fut/cgm_test \
	lib/github.com/jonesz/krylov-fut/*.c \
	lib/github.com/jonesz/krylov-fut/*.fut.* \

.PHONY: check test clean

all: check test

check:
	futhark check lib/github.com/jonesz/krylov-fut/cgm.fut

test:
	futhark test lib/github.com/jonesz/krylov-fut/cgm_test.fut

clean:	
	$(RM) $(TEST_FRAGMENTS)
