import argparse
import sys

def my_func_that_return_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('foo', default=False, help='foo help')
	parser.add_argument('bar', default=False)

	subparsers = parser.add_subparsers()

	subparser = subparsers.add_parser('install', help='install help')
	subparser.add_argument('ref', type=str, help='foo1 help')
	subparser.add_argument('--upgrade', action='store_true', default=False, help='foo2 help')

#	if len(sys.argv) == 1:
#		print >> sys.stderr,parser.print_help()
#		exit(0)
#	args = parser.parse_args()
	return parser
