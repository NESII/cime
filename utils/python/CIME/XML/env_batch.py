"""
Interface to the env_batch.xml file.  This class inherits from EnvBase
"""
import stat
import time
from CIME.XML.standard_module_setup import *
from CIME.task_maker import TaskMaker
from CIME.utils import convert_to_type
from CIME.XML.env_base import EnvBase
from CIME.utils import convert_to_string, transform_vars, get_cime_root
from copy import deepcopy

import re

logger = logging.getLogger(__name__)

class EnvBatch(EnvBase):

    def __init__(self, case_root=None, infile="env_batch.xml"):
        """
        initialize an object interface to file env_batch.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
        self.prereq_jobid = None
        self.batchtype = None

    def set_value(self, item, value, subgroup=None, ignore_type=False):
        val = None
        if item == "JOB_WALLCLOCK_TIME":
            # Most systems use %H:%M:%S format for wallclock but LSF
            # uses %H:%M this code corrects the value passed in to be
            # the correct format - if we find we have more exceptions
            # than this we may need to generalize this further
            walltime_format = self.get_value("walltime_format", subgroup=None)
            if walltime_format is not None and walltime_format.count(":") != value.count(":"):
                t = time.strptime(value,"%H:%M:%S")
                value = time.strftime(walltime_format, t)

        # allow the user to set all instances of item if subgroup is not provided
        if subgroup is None:
            nodes = self.get_nodes("entry", {"id":item})
            for node in nodes:
                self._set_value(node, value, vid=item, ignore_type=ignore_type)
                val = value
        else:
            nodes = self.get_nodes("job", {"name":subgroup})
            for node in nodes:
                vnode = self.get_optional_node("entry", {"id":item}, root=node)
                if vnode is not None:
                    val = self._set_value(vnode, value, vid=item, ignore_type=ignore_type)

        return val

    def get_value(self, item, attribute=None, resolved=True, subgroup="case.run"):
        """
        Must default subgroup to something in order to provide single return value
        """
        value = None
        if subgroup is None:
            node = self.get_optional_node(item, attribute)
            if node is not None:
                value = node.text
                if resolved:
                    value = self.get_resolved_value(value)
            else:
                value = EnvBase.get_value(self,item,attribute,resolved)
        else:
            job_node = self.get_optional_node("job", {"name":subgroup})
            if job_node is not None:
                node = self.get_optional_node("entry", {"id":item}, root=job_node)
                if node is not None:
                    value = self.get_resolved_value(node.get("value"))

                    # Return value as right type if we were able to fully resolve
                    # otherwise, we have to leave as string.
                    if "$" not in value:
                        type_str = self._get_type_info(node)
                        value = convert_to_type(value, type_str, item)
        return value

    def get_values(self, item, attribute=None, resolved=True, subgroup=None):
        """Returns the value as a string of the first xml element with item as attribute value.
        <elememt_name attribute='attribute_value>value</element_name>"""

        logger.debug("(get_values) Input values: %s , %s , %s , %s , %s" , self.__class__.__name__ , item, attribute, resolved, subgroup)

        nodes   = [] # List of identified xml elements
        results = [] # List of identified parameters
             
        
        # Find all nodes with attribute name and attribute value item
        # xpath .//*[name='item']
        # for job in self.get_nodes("job") :
        
        groups = self.get_nodes("group")

        for group in groups :

            roots = []
            jobs  = []
            jobs  = self.get_nodes("job" , root=group)

            if (len(jobs)) :
                roots = jobs
            else :
                roots = [group]

            for root in roots :

                if item :
                    nodes = self.get_nodes("entry",{"id" : item} , root=root )
                else :
                    # Return all nodes
                    nodes = self.get_nodes("entry" , root=root)

                # seach in all entry nodes
                for node in nodes:

                    
                    # Build return structure
                    attr          = node.get('id')
                    group_name     = None
                    
                    # determine group
                    if (root.tag == "job") :
                        group_name = root.get('name') 
                    else:
                        group_name = root.get('id')

                    val             = node.get('value')
                    attribute_type  = self._get_type(node)
                    desc            = self._get_description(node)
                    default         = super(EnvBatch , self)._get_default(node)
                    filename        = self.filename

                    tmp = { 'group' : group_name , 'attribute' : attr , 'value' : val , 'type' : attribute_type , 'description' : desc , 'default' : default , 'file' : filename}
                    logger.debug("Found node with value for %s = %s" , item , tmp )
                    
                    # add single result to list
                    results.append(tmp) 

        logger.debug("(get_values) Return value:  %s" , results )

        return results

    def get_type_info(self, vid):
        nodes = self.get_nodes("entry",{"id":vid})
        type_info = None
        for node in nodes:
            new_type_info = self._get_type_info(node)
            if type_info is None:
                type_info = new_type_info
            else:
                expect( type_info == new_type_info,
                        "Inconsistent type_info for entry id=%s %s %s" % (vid, new_type_info, type_info))
        return type_info

    def get_jobs(self):
        jobs = []
        for node in self.get_nodes("job"):
            name = node.get("name")
            jobs.append(name)
        return jobs

    def create_job_groups(self, bjobs):
        # only the job_submission group is repeated
        group = self.get_node("group", {"id":"job_submission"})
        # look to see if any jobs are already defined
        cjobs = self.get_jobs()
        childnodes = []

        expect(len(cjobs)==0," Looks like job groups have already been created")

        for child in reversed(group):
            childnodes.append(deepcopy(child))
            group.remove(child)

        for name,jdict in bjobs:
            newjob = ET.Element("job")
            newjob.set("name",name)
            for field in jdict.keys():
                val = jdict[field]
                node = ET.SubElement(newjob, "entry", {"id":field,"value":val})
                tnode = ET.SubElement(node, "type")
                tnode.text = "char"
            for child in childnodes:
                newjob.append(deepcopy(child))
            group.append(newjob)

    def cleanupnode(self, node):
        if node.get("id") == "batch_system":
            fnode = node.find(".//file")
            node.remove(fnode)
            gnode = node.find(".//group")
            node.remove(gnode)
            vnode = node.find(".//values")
            if vnode is not None:
                node.remove(vnode)
        else:
            node = EnvBase.cleanupnode(self, node)
        return node

    def set_batch_system(self, batchobj, batch_system_type=None):
        if batch_system_type is not None:
            self.set_batch_system_type(batch_system_type)
        if batchobj.batch_system_node is not None:
            self.root.append(deepcopy(batchobj.batch_system_node))
        if batchobj.machine_node is not None:
            self.root.append(deepcopy(batchobj.machine_node))

    def make_batch_script(self, input_template, job, case):
        expect(os.path.exists(input_template), "input file '%s' does not exist" % input_template)

        task_maker = TaskMaker(case)

        self.maxthreads = task_maker.maxthreads
        self.taskgeometry = task_maker.taskgeometry
        self.threadgeometry = task_maker.threadgeometry
        self.taskcount = task_maker.taskcount
        self.thread_count = task_maker.thread_count
        self.pedocumentation = task_maker.document()
        self.ptile = task_maker.ptile
        self.tasks_per_node = task_maker.tasks_per_node
        self.max_tasks_per_node = task_maker.MAX_TASKS_PER_NODE
        self.tasks_per_numa = task_maker.tasks_per_numa
        self.num_tasks = task_maker.totaltasks

        task_count = self.get_value("task_count")
        if task_count == "default":
            self.sumpes = task_maker.fullsum
            self.totaltasks = task_maker.totaltasks
            self.fullsum = task_maker.fullsum
            self.sumtasks = task_maker.totaltasks
            self.task_count = task_maker.fullsum
            self.num_nodes = task_maker.num_nodes
        else:
            self.sumpes = task_count
            self.totaltasks = task_count
            self.fullsum = task_count
            self.sumtasks = task_count
            self.task_count = task_count
            self.num_nodes = task_count
            self.pedocumentation = ""
        self.job_id = case.get_value("CASE") + os.path.splitext(job)[1]
        if "pleiades" in case.get_value("MACH"):
            # pleiades jobname needs to be limited to 15 chars
            self.job_id = self.job_id[:15]
        self.output_error_path = self.job_id

        self.batchdirectives = self.get_batch_directives(case, job)

        output_text = transform_vars(open(input_template,"r").read(), case=case, subgroup=job, check_members=self)
        with open(job, "w") as fd:
            fd.write(output_text)
        os.chmod(job, os.stat(job).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def set_job_defaults(self, bjobs, pesize=None, walltime=None):
        if self.batchtype is None:
            self.batchtype = self.get_batch_system_type()
        if self.batchtype == 'none':
            return
        for job, jsect in bjobs:
            task_count = jsect["task_count"]
            if task_count is None or task_count == "default":
                task_count = pesize
            else:
                task_count = int(task_count)
            queue = self.select_best_queue(task_count)
            self.set_value("JOB_QUEUE", queue, subgroup=job)
            walltime = self.get_max_walltime(queue) if walltime is None else walltime
            if walltime is None:
                walltime = self.get_default_walltime()
            self.set_value( "JOB_WALLCLOCK_TIME", walltime , subgroup=job)
            logger.info("Job %s queue %s walltime %s"%(job, queue, walltime))

    def get_batch_directives(self, case, job, raw=False):
        """
        """
        result = []
        directive_prefix = self.get_node("batch_directive").text
        directive_prefix = "" if directive_prefix is None else directive_prefix

        roots = self.get_nodes("batch_system")
        for root in roots:
            if root is not None:
                nodes = self.get_nodes("directive", root=root)
                for node in nodes:
                    directive = self.get_resolved_value("" if node.text is None else node.text)
                    default = node.get("default")
                    if not raw:
                        directive = transform_vars(directive, case=case, subgroup=job, default=default, check_members=self)
                    result.append("%s %s" % (directive_prefix, directive))

        return "\n".join(result)

    def get_submit_args(self, case, job):
        '''
        return a list of touples (flag, name)
        '''
        submitargs = " "
        bs_nodes = self.get_nodes("batch_system")
        submit_arg_nodes = []
        for node in bs_nodes:
            submit_arg_nodes += self.get_nodes("arg",root=node)
        for arg in submit_arg_nodes:
            flag = arg.get("flag")
            name = arg.get("name")
            if name is None:
                submitargs+=" %s"%flag
            else:
                val = case.get_value(name,subgroup=job)
                if val is None:
                    val = case.get_resolved_value(name)

                if val is not None and len(val) > 0 and val != "None":
                    if flag.rfind("=", len(flag)-1, len(flag)) >= 0 or\
                       flag.rfind(":", len(flag)-1, len(flag)) >= 0:
                        submitargs+=" %s%s"%(flag,str(val).strip())
                    else:
                        submitargs+=" %s %s"%(flag,str(val).strip())

        return submitargs

    def submit_jobs(self, case, no_batch=False, job=None):
        alljobs = self.get_jobs()
        startindex = 0
        jobs = []
        if job is not None:
            expect(job in alljobs, "Do not know about batch job %s"%job)
            startindex = alljobs.index(job)

        for index, job in enumerate(alljobs):
            if index < startindex:
                continue
            logger.debug( "Index %d job %s"%(index, job))
            try:
                prereq = case.get_resolved_value(self.get_value('prereq', subgroup=job))
                prereq = eval(prereq)
            except:
                expect(False,"Unable to evaluate prereq expression '%s' for job '%s'"%(self.get_value('prereq',subgroup=job), job))
            if prereq:
                jobs.append((job,self.get_value('dependency', subgroup=job)))
        depid = {}
        for job, dependency in jobs:
            if dependency is not None:
                deps = dependency.split()
            else:
                deps = []
            jobid = ""
            if self.prereq_jobid is not None:
                jobid = self.prereq_jobid
            for dep in deps:
                if dep in depid.keys() and depid[dep] is not None:
                    jobid += " "+str(depid[dep])
#TODO: doubt these will be used
#               elif dep == "and":
#                   jobid += " && "
#               elif dep == "or":
#                   jobid += " || "


            slen = len(jobid)
            if slen == 0:
                jobid = None

            depid[job] = self.submit_single_job(case, job, jobid, no_batch=no_batch)

    def submit_single_job(self, case, job, depid=None, no_batch=False):
        caseroot = case.get_value("CASEROOT")
        batch_system = self.get_value("BATCH_SYSTEM", subgroup=None)
        if batch_system is None or batch_system == "none" or no_batch:
            # Import here to avoid circular include
            from CIME.case_test       import case_test
            from CIME.case_run        import case_run
            from CIME.case_st_archive import case_st_archive
            from CIME.case_lt_archive import case_lt_archive

            logger.info("Starting job script %s" % job)

            # Hack until all testcases are ported to python
            testcase = case.get_value("TESTCASE")
            cimeroot = get_cime_root()
            testscript = os.path.join(cimeroot, "scripts", "Testing", "Testcases", "%s_script" % testcase)
            if job == "case.test" and testcase is not None and os.path.exists(testscript):
                run_cmd("%s --caseroot %s" % (os.path.join(".", job), caseroot))
            else:
                # This is what we want longterm
                function_name = job.replace(".", "_")
                locals()[function_name](case)

            return

        submitargs = self.get_submit_args(case, job)

        if depid is not None:
            dep_string = self.get_value("depend_string", subgroup=None)
            dep_string = dep_string.replace("jobid",depid.strip())
            submitargs += " " + dep_string

        batchsubmit = self.get_value("batch_submit", subgroup=None)
        expect(batchsubmit is not None,
               "Unable to determine the correct command for batch submission.")
        batchredirect = self.get_value("batch_redirect", subgroup=None)
        submitcmd = ''
        for string in (batchsubmit, submitargs, batchredirect, job):
            if  string is not None:
                submitcmd += string + " "

        if self.batchtype == "pbs":
            submitcmd += " -F \"--caseroot %s\""%caseroot

        logger.info("Submitting job script %s"%submitcmd)
        output = run_cmd(submitcmd)
        jobid = self.get_job_id(output)
        logger.debug("Submitted job id is %s"%jobid)
        return jobid

    def get_batch_system_type(self):
        nodes = self.get_nodes("batch_system")
        for node in nodes:
            type_ = node.get("type")
            if type_ is not None:
                self.batchtype = type_
        return self.batchtype

    def set_batch_system_type(self, batchtype):
        self.batchtype = batchtype

    def get_job_id(self, output):
        jobid_pattern = self.get_value("jobid_pattern", subgroup=None)
        expect(jobid_pattern is not None, "Could not find jobid_pattern in env_batch.xml")
        jobid = re.search(jobid_pattern, output).group(1)
        return jobid

    def select_best_queue(self, num_pes):
        # Make sure to check default queue first.
        all_queues = []
        all_queues.append( self.get_default_queue())
        all_queues = all_queues + self.get_all_queues()
        for queue in all_queues:
            if queue is not None:
                jobmin = queue.get("jobmin")
                jobmax = queue.get("jobmax")
                # if the fullsum is between the min and max # jobs, then use this queue.
                if jobmin is not None and jobmax is not None and num_pes >= int(jobmin) and num_pes <= int(jobmax):
                    return queue.text
        return None

    def get_max_walltime(self, queue):
        for queue_node in self.get_all_queues():
            if queue_node.text == queue:
                return queue_node.get("walltimemax")

    def get_default_walltime(self):
        return self.get_value("walltime", attribute={"default" : "true"}, subgroup=None)

    def get_default_queue(self):
        return self.get_optional_node("queue", attributes={"default" : "true"})

    def get_all_queues(self):
        return self.get_nodes("queue")