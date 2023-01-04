TEST_FILES = \
	lib/github.com/jonesz/krylov-fut/symmetric_test.fut \
	lib/github.com/jonesz/krylov-fut/cgm_test.fut

TEST_JUNK = \
	lib/github.com/jonesz/krylov-fut/symmetric_test \
	lib/github.com/jonesz/krylov-fut/cgm_test \
	lib/github.com/jonesz/krylov-fut/*.c \
	lib/github.com/jonesz/krylov-fut/*.fut.* \

.PHONY: clean test

test:
	futhark test $(TEST_FILES)

clean:	
	$(RM) $(TEST_JUNK)
