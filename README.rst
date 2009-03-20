A regression test suite for OpenBabel
=====================================

Yes, OpenBabel already has a test suite, but the more tests the better.

This test suite is written in Python, and currently focuses on regressions (although any test is welcome). The current tests are named after the entries in OpenBabel's bug tracker (go to http://openbabel.sf.net, click on Trackers, then Bugs).

How to download
---------------

One way is to click the download button at http://github.com/baoilleach/obunittest/tree/master.

A better way is to install Git. On Windows, you can install Git under Cygwin. Then use the following command::

   git clone git://github.com/baoilleach/obunittest.git

With Git, you can stay up to the date with the latest version by just typing "git fetch".

How to run
----------

To run all tests::

   python sweet.py

To run one test::

   python sweet.py PR1733905.py

How to add a test
-----------------

Copy an old test, rename it, edit it, and email me the result.

A better way is to fork the repository on Github, "git clone" your fork, make your changes, "git push" it to github, and then send me a pull request.

What's a good test?
-------------------

A poor test is very specific to the bug or feature being tested. It tests nothing else and will only fail if a very specific problem occurs. 

A good test maximises the chances of failing. It should fail if the specific problem occurs, but it should also fail if any 1 of a 100 things goes wrong.

From the point of view of maintaining tests, a good test should not rely on anything that may change in the future. If a test fails, it should indicate an error in OpenBabel, not an error in the test.

How to serialise a molecule
---------------------------

Given a molecule in a file, you can create a serialised version using::

   python sweet.py serial myfile.mol

Copy and paste the result into your test file to ensure that the molecule's atom types and bond orders are correctly perceived in all future OpenBabel versions. See testIdentity in PR1739905.py for an example of usage.
